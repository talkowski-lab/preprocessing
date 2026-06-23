version 1.0

## Portions Copyright Broad Institute, 2018
##
## This WDL pipeline implements QC in human whole-genome or exome/targeted sequencing data.
##
## Requirements/expectations
## - Human paired-end sequencing data in aligned BAM or CRAM format
## - Input BAM/CRAM files must additionally comply with the following requirements:
## - - files must pass validation by ValidateSamFile
## - - reads are provided in query-sorted order
## - - all reads must have an RG tag
## - Reference genome must be Hg38 with ALT contigs
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL open source code license (BSD-3).
## Full license text at https://github.com/openwdl/wdl/blob/master/LICENSE
## Note however that the programs it calls may be subject to different licenses. 
## Users are responsible for checking that they are authorized to run all programs before running this script.
## - [Picard](https://broadinstitute.github.io/picard/)
## - [VerifyBamID2](https://github.com/Griffan/VerifyBamID)

# Git URL import
#import "tasks/Qc.wdl" as QC

# WORKFLOW DEFINITION
workflow SingleSampleQc {
  input {
    File input_bam
    File input_bam_index
    File input_gvcf
    File input_gvcf_index
    File dbsnp_vcf
    File dbsnp_vcf_index
    File ref_cache
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String base_name
    Int preemptible_tries
    File coverage_interval_list
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File evaluation_interval_list
    Boolean is_wgs
    Boolean? is_outlier_data
    Boolean is_gvcf = true

    File evaluation_thresholds
  }

  # Not overridable:
  Int read_length = 250
  
  # Generate a BAM or CRAM index
  # Estimate level of cross-sample contamination
  call CheckContamination {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      contamination_sites_ud = contamination_sites_ud,
      contamination_sites_bed = contamination_sites_bed,
      contamination_sites_mu = contamination_sites_mu,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      output_prefix = base_name + ".verify_bam_id",
      preemptible_tries = preemptible_tries,
  }

  # Calculate the duplication rate since MarkDuplicates was already performed
  call CollectDuplicateMetrics {
    input:
      input_bam = input_bam,
      input_bam_index = input_bam_index,
      output_bam_prefix = base_name,
      ref_dict = ref_dict,
      ref_fasta = ref_fasta,
      ref_fasta_index = ref_fasta_index,
      preemptible_tries = preemptible_tries
  }

call CollectVariantCallingMetrics {
    input:
      input_vcf = input_gvcf,
      input_vcf_index = input_gvcf_index,
      metrics_basename = base_name,
      dbsnp_vcf = dbsnp_vcf,
      dbsnp_vcf_index = dbsnp_vcf_index,
      ref_dict = ref_dict,
      evaluation_interval_list = evaluation_interval_list,
      is_gvcf = is_gvcf,
      preemptible_tries = preemptible_tries
  }

  # Outputs that will be retained when execution is complete
  output {


    File selfSM = CheckContamination.selfSM
    Float contamination = CheckContamination.contamination

    File duplication_metrics_file = CollectDuplicateMetrics.duplication_metrics_file
    String percent_duplication = CollectDuplicateMetrics.percent_duplication

    File vcf_summary_metrics = CollectVariantCallingMetrics.summary_metrics
    File vcf_detail_metrics = CollectVariantCallingMetrics.detail_metrics

  }
}
task CheckContamination {
  input {
    File input_bam
    File input_bam_index
    File contamination_sites_ud
    File contamination_sites_bed
    File contamination_sites_mu
    File ref_fasta
    File ref_fasta_index
    String output_prefix
    Int preemptible_tries
    Boolean disable_sanity_check = false
  }

  Int disk_size = ceil(size(input_bam, "GiB") + size(ref_fasta, "GiB")) + 30

  command <<<
    set -e

    # creates a ~{output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    /usr/gitc/VerifyBamID \
    --Verbose \
    --NumPC 4 \
    --Output ~{output_prefix} \
    --BamFile ~{input_bam} \
    --Reference ~{ref_fasta} \
    --UDPath ~{contamination_sites_ud} \
    --MeanPath ~{contamination_sites_mu} \
    --BedPath ~{contamination_sites_bed} \
    ~{true="--DisableSanityCheck" false="" disable_sanity_check} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('~{output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"]))
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    memory: "4 GiB"
    disks: "local-disk " + disk_size + " HDD"
    docker: "us.gcr.io/broad-gotc-prod/verify-bam-id:c1cba76e979904eb69c31520a0d7f5be63c72253-1553018888"
    cpu: "2"
  }
  output {
    File selfSM = "~{output_prefix}.selfSM"
    Float contamination = read_float(stdout())
    Map[String, String] metrics = { "FREEMIX": read_string(stdout()) }
  }
}

task CollectVariantCallingMetrics {
  input {
    File input_vcf
    File input_vcf_index
    String metrics_basename
    File dbsnp_vcf
    File dbsnp_vcf_index
    File ref_dict
    File evaluation_interval_list
    Boolean is_gvcf = true
    Int preemptible_tries
  }

  Int disk_size = ceil(size(input_vcf, "GiB") + size(dbsnp_vcf, "GiB")) + 20

  command {
    java -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
      CollectVariantCallingMetrics \
      INPUT=~{input_vcf} \
      OUTPUT=~{metrics_basename} \
      DBSNP=~{dbsnp_vcf} \
      SEQUENCE_DICTIONARY=~{ref_dict} \
      TARGET_INTERVALS=~{evaluation_interval_list} \
      ~{true="GVCF_INPUT=true" false="" is_gvcf}
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    preemptible: preemptible_tries
    memory: "3000 MiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File summary_metrics = "~{metrics_basename}.variant_calling_summary_metrics"
    File detail_metrics = "~{metrics_basename}.variant_calling_detail_metrics"
  }
}

task CollectDuplicateMetrics {
  input {
    File input_bam
    File input_bam_index
    String output_bam_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int preemptible_tries
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil(size(input_bam, "GiB") + ref_size) + 20

  String duplication_metric_object_file = '~{output_bam_prefix}.duplication_metrics.metrics_only'

  command <<<
    java -Xms5000m -jar /usr/picard/picard.jar \
      CollectDuplicateMetrics \
      METRICS_FILE=~{output_bam_prefix}.duplication_metrics \
      INPUT=~{input_bam} \
      ASSUME_SORTED=true \
      REFERENCE_SEQUENCE=~{ref_fasta}

    grep -v '#' '~{output_bam_prefix}.duplication_metrics' | grep '.\+' | perl -E 'my ($keys, $values) = <>; chomp $keys; chomp $values; my @k = split("\t", $keys); my @v = split("\t", $values); for(0..$#k) { say join("\t", $k[$_], $v[$_]); }'      > ~{duplication_metric_object_file}

  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/picard-cloud:2.21.7"
    memory: "7 GiB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: preemptible_tries
  }
  output {
    File duplication_metrics_file = "~{output_bam_prefix}.duplication_metrics"
    Map[String, String] duplication_metrics = read_map(duplication_metric_object_file)
    String percent_duplication = duplication_metrics["PERCENT_DUPLICATION"]
  }
}


