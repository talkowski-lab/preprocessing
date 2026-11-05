version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow MergeVCFSamplesUpdateInfo {
    input {
        Array[File] vcf_files
        String merge_vcf_samples_update_info_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/main/scripts/hail_merge_vcf_samples_update_info.py"
        String hail_vcf_utils_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/tags/hail-vcf-utils-v0.1/scripts/hail_vcf_utils.py"
        String hail_docker
        String prefix
        String genome_build = "GRCh38"

        Int batch_size = 10
        String? output_vcf_filename_override  # End with ".bgz" or ".gz" (output is bgzipped and indexed)
        Array[String] min_info_fields = ["FractionInformativeReads", "LOD", "MQ", "MQRankSum", "BaseQRankSum", "ReadPosRankSum", "QD", "VQSLOD", "QUALapprox", "ExcessHet"]
        Array[String] max_info_fields = ["FS", "SOR"]
        Array[String] sum_info_fields = ["DP", "MQ_DP", "VarDP"]

        # false: process each VCF on its own (subset/rename/FORMAT fields/AC-AF-AN), without merging
        Boolean merge_inputs = true
        # OPTIONAL: one sample ID per line; only these samples are kept
        File? sample_file
        # OPTIONAL: one prefix per VCF, in vcf_files order; samples are renamed to <prefix>:<sample>
        Array[String] cohort_prefixes = []
        # OPTIONAL: FORMAT fields to drop (GT cannot be dropped); empty keeps every field shared by all VCFs
        Array[String] format_fields_to_drop = []
        # outer: keep variants in any VCF; inner: keep only variants in every VCF
        String site_join = "outer"
        # genotypes for samples at a variant only in other VCFs:
        # all: every missing GT becomes 0/0, including no-calls (./.) in the input VCFs
        # absent: 0/0 only where the sample's VCF has no record for the variant
        # absent_at_genotyped_loci: as absent, but only where the sample's VCF has a record at the same locus (e.g. another allele)
        # none: missing GTs stay missing
        # GATK already writes samples with no reads as 0/0 (GQ 0, DP 0), so with all/absent AN counts them too;
        # with absent_at_genotyped_loci/none, a variant in only one VCF gets AN from that VCF's samples alone, inflating AF.
        # with all/absent, recompute AN/AC/AF downstream from genotypes that pass QC (e.g. GQ, DP, AB cutoffs)
        # every call is made diploid before merging (haploid DRAGEN male chrX/chrY calls become 0/0 or 1/1), so fills are 0/0
        String fill_missing = "absent"
        # OPTIONAL: VCF header to use for every input VCF
        File? header_file

        RuntimeAttr? runtime_attr_merge
    }

    if (merge_inputs) {
        String output_vcf_file = select_first([output_vcf_filename_override, prefix + '.merged.vcf.bgz'])

        call mergeVCFSamplesUpdateInfo {
            input:
                merge_vcf_samples_update_info_script = merge_vcf_samples_update_info_script,
                hail_vcf_utils_script = hail_vcf_utils_script,
                vcf_files = vcf_files,
                batch_size = batch_size,
                output_vcf_file = output_vcf_file,
                min_info_fields = min_info_fields,
                max_info_fields = max_info_fields,
                sum_info_fields = sum_info_fields,
                genome_build = genome_build,
                sample_file = sample_file,
                cohort_prefixes = cohort_prefixes,
                format_fields_to_drop = format_fields_to_drop,
                site_join = site_join,
                fill_missing = fill_missing,
                header_file = header_file,
                hail_docker = hail_docker,
                runtime_attr_override = runtime_attr_merge
        }
    }

    if (!merge_inputs) {
        scatter (i in range(length(vcf_files))) {
            call mergeVCFSamplesUpdateInfo as processVCF {
                input:
                    merge_vcf_samples_update_info_script = merge_vcf_samples_update_info_script,
                    hail_vcf_utils_script = hail_vcf_utils_script,
                    vcf_files = [vcf_files[i]],
                    batch_size = batch_size,
                    output_vcf_file = sub(basename(vcf_files[i]), "\\.vcf(\\.b?gz)?$", "") + '.processed.vcf.bgz',
                    min_info_fields = min_info_fields,
                    max_info_fields = max_info_fields,
                    sum_info_fields = sum_info_fields,
                    genome_build = genome_build,
                    sample_file = sample_file,
                    cohort_prefixes = if length(cohort_prefixes) > 0 then [cohort_prefixes[i]] else [],
                    format_fields_to_drop = format_fields_to_drop,
                    site_join = site_join,
                    fill_missing = fill_missing,
                    header_file = header_file,
                    hail_docker = hail_docker,
                    runtime_attr_override = runtime_attr_merge
            }
        }
    }

    output {
        File? merged_vcf_file = mergeVCFSamplesUpdateInfo.merged_vcf_file
        File? merged_vcf_idx = mergeVCFSamplesUpdateInfo.merged_vcf_idx
        Array[File]? processed_vcf_files = processVCF.merged_vcf_file
        Array[File]? processed_vcf_idx = processVCF.merged_vcf_idx
    }
}

task mergeVCFSamplesUpdateInfo {
    input {
        String merge_vcf_samples_update_info_script
        String hail_vcf_utils_script
        Array[File] vcf_files
        Int batch_size
        String output_vcf_file
        String hail_docker
        Array[String] min_info_fields
        Array[String] max_info_fields
        Array[String] sum_info_fields
        String genome_build
        File? sample_file
        Array[String] cohort_prefixes
        Array[String] format_fields_to_drop
        String site_join
        String fill_missing
        File? header_file
        Float input_disk_scale = 5.0
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, 'GB')
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
        mem_gb: 4,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])

    runtime {
        memory: "~{memory} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: cpu_cores
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -e

        # Write the Python script to a temporary file
        curl -L ~{merge_vcf_samples_update_info_script} -o merge_vcf_files.py
        curl -L ~{hail_vcf_utils_script} -o hail_vcf_utils.py

        # Run the Python script with the provided arguments
        python3 merge_vcf_files.py \
            --vcf-files ~{sep=" " vcf_files} \
            --batch-size ~{batch_size} \
            --output-vcf-file ~{output_vcf_file} \
            --mem ~{memory} \
            --min-info-fields ~{sep=" " min_info_fields} \
            --max-info-fields ~{sep=" " max_info_fields} \
            --sum-info-fields ~{sep=" " sum_info_fields} \
            --genome-build ~{genome_build} \
            --cohort-prefixes ~{sep=" " cohort_prefixes} \
            --format-fields-to-drop ~{sep=" " format_fields_to_drop} \
            --site-join ~{site_join} \
            --fill-missing ~{fill_missing} \
            ~{"--sample-file " + sample_file} \
            ~{"--header-file " + header_file}
    >>>

    output {
        File merged_vcf_file = output_vcf_file
        File merged_vcf_idx = output_vcf_file + ".tbi"
    }
}
