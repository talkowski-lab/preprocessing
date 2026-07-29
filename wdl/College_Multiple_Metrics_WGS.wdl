workflow CollectMultipleMetricsWgs {
    File input_bam
    File input_bam_index
    String base_name = basename(input_bam, ".bam")
    
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File wgs_coverage_interval_list

    # average read length in the file: Picard's default 150bp
    Int? read_length_input
    Int read_length = select_first([read_length_input, 150])

    # optional inputs for runtime parameters
    Int? preemptible_attempts
    Int? max_retries
    Int? disk_pad

    # fixed runtime parameters
    String? docker_override
    String docker = select_first([docker_override, "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.2-1510681135"])
    Int preemptible_tries = select_first([preemptible_attempts, 3])
    Int max_retry_param = select_first([max_retries, 2])
    Int disk_size = ceil(bam_size + ref_size) + select_first([disk_pad, 10])

    # disk size parameters
    Float bam_size = size(input_bam, "GB") + size(input_bam_index, "GB")
    Float ref_size = size(ref_fasta, "GB") + size(ref_fasta_index, "GB") + size(ref_dict, "GB")

    call CollectReadgroupBamQualityMetrics {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        output_prefix = base_name + ".readgroup",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = docker,
        preemptible_tries = preemptible_tries,
        max_retries = max_retry_param,
        disk_size = disk_size
    }

    call CollectAggregationMetrics {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        output_prefix = base_name,
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = docker,
        preemptible_tries = preemptible_tries,
        max_retries = max_retry_param,
        disk_size = disk_size
    }

    call CollectWgsMetrics {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        output_prefix = base_name + ".wgs_metrics",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        wgs_coverage_interval_list = wgs_coverage_interval_list,
        read_length = read_length,
        docker = docker,
        preemptible_tries = preemptible_tries,
        max_retries = max_retry_param,
        disk_size = disk_size
    }

    call CollectRawWgsMetrics {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        output_prefix = base_name + ".raw_wgs_metrics",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        wgs_coverage_interval_list = wgs_coverage_interval_list,
        read_length = read_length,
        docker = docker,
        preemptible_tries = preemptible_tries,
        max_retries = max_retry_param,
        disk_size = disk_size
    }

    call ValidateSamFile {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        output_prefix = base_name + ".validation_report",
        ref_dict = ref_dict,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        docker = docker,
        preemptible_tries = preemptible_tries,
        max_retries = max_retry_param,
        disk_size = disk_size
    }

    output {
        File readgroup_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics
        File readgroup_gc_bias_detail_metrics = CollectReadgroupBamQualityMetrics.gc_bias_detail_metrics
        File readgroup_gc_bias_pdf = CollectReadgroupBamQualityMetrics.gc_bias_pdf
        File readgroup_gc_bias_summary_metrics = CollectReadgroupBamQualityMetrics.gc_bias_summary_metrics

        File alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics
        File bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics
        File bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics
        File gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
        File gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
        File gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics
        File insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf
        File insert_size_metrics = CollectAggregationMetrics.insert_size_metrics
        File pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics
        File pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics
        File quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
        File quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics

        File wgs_metrics = CollectWgsMetrics.metrics
        File raw_wgs_metrics = CollectRawWgsMetrics.metrics
        File validate_sam_report = ValidateSamFile.report
    }
}



# Collect alignment summary and GC bias quality metrics
task CollectReadgroupBamQualityMetrics {
    File input_bam
    File input_bam_index
    String output_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String docker
    Int preemptible_tries
    Int max_retries
    Int disk_size

    command {
        java -Xms5000m -jar /usr/gitc/picard.jar \
            CollectMultipleMetrics \
            INPUT=${input_bam} \
            REFERENCE_SEQUENCE=${ref_fasta} \
            OUTPUT=${output_prefix} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            PROGRAM="CollectGcBiasMetrics" \
            METRIC_ACCUMULATION_LEVEL="null" \
            METRIC_ACCUMULATION_LEVEL="READ_GROUP"
    }

    runtime {
        docker: docker
        memory: "7 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_tries
        maxRetries: max_retries
    }
    output {
        File alignment_summary_metrics = "${output_prefix}.alignment_summary_metrics"
        File gc_bias_detail_metrics = "${output_prefix}.gc_bias.detail_metrics"
        File gc_bias_pdf = "${output_prefix}.gc_bias.pdf"
        File gc_bias_summary_metrics = "${output_prefix}.gc_bias.summary_metrics"
    }
}

# Collect quality metrics from the aggregated bam
task CollectAggregationMetrics {
    File input_bam
    File input_bam_index
    String output_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String docker
    Int preemptible_tries
    Int max_retries
    Int disk_size

    command {
        java -Xms5000m -jar /usr/gitc/picard.jar \
            CollectMultipleMetrics \
            INPUT=${input_bam} \
            REFERENCE_SEQUENCE=${ref_fasta} \
            OUTPUT=${output_prefix} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            PROGRAM="CollectInsertSizeMetrics" \
            PROGRAM="CollectSequencingArtifactMetrics" \
            PROGRAM="CollectGcBiasMetrics" \
            PROGRAM="QualityScoreDistribution" \
            METRIC_ACCUMULATION_LEVEL="null" \
            METRIC_ACCUMULATION_LEVEL="SAMPLE" \
            METRIC_ACCUMULATION_LEVEL="LIBRARY"

        touch ${output_prefix}.insert_size_metrics
        touch ${output_prefix}.insert_size_histogram.pdf
    }
    runtime {
        docker: docker
        memory: "7 GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: preemptible_tries
        maxRetries: max_retries
    }
    output {
        File alignment_summary_metrics = "${output_prefix}.alignment_summary_metrics"
        File bait_bias_detail_metrics = "${output_prefix}.bait_bias_detail_metrics"
        File bait_bias_summary_metrics = "${output_prefix}.bait_bias_summary_metrics"
        File gc_bias_detail_metrics = "${output_prefix}.gc_bias.detail_metrics"
        File gc_bias_pdf = "${output_prefix}.gc_bias.pdf"
        File gc_bias_summary_metrics = "${output_prefix}.gc_bias.summary_metrics"
        File insert_size_histogram_pdf = "${output_prefix}.insert_size_histogram.pdf"
        File insert_size_metrics = "${output_prefix}.insert_size_metrics"
        File pre_adapter_detail_metrics = "${output_prefix}.pre_adapter_detail_metrics"
        File pre_adapter_summary_metrics = "${output_prefix}.pre_adapter_summary_metrics"
        File quality_distribution_pdf = "${output_prefix}.quality_distribution.pdf"
        File quality_distribution_metrics = "${output_prefix}.quality_distribution_metrics"
    }
}

# Note these tasks will break if the read lengths in the bam are greater than 250.
task CollectWgsMetrics {
    File input_bam
    File input_bam_index
    String output_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    
    File wgs_coverage_interval_list
    Int read_length

    String docker
    Int preemptible_tries
    Int max_retries
    Int disk_size

    command {
        java -Xms2000m -jar /usr/gitc/picard.jar \
            CollectWgsMetrics \
            INPUT=${input_bam} \
            VALIDATION_STRINGENCY=SILENT \
            REFERENCE_SEQUENCE=${ref_fasta} \
            INCLUDE_BQ_HISTOGRAM=true \
            INTERVALS=${wgs_coverage_interval_list} \
            OUTPUT=${output_prefix} \
            USE_FAST_ALGORITHM=true \
            READ_LENGTH=${read_length}
    }
    runtime {
        docker: docker
        preemptible: preemptible_tries
        memory: "3 GB"
        disks: "local-disk " + disk_size + " HDD"
        maxRetries: max_retries
    }
    output {
        File metrics = "${output_prefix}"
    }
}

# Collect raw WGS metrics (commonly used QC thresholds)
task CollectRawWgsMetrics {
    File input_bam
    File input_bam_index
    String output_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    
    File wgs_coverage_interval_list
    Int read_length

    String docker
    Int preemptible_tries
    Int max_retries
    Int disk_size

    command {
        java -Xms2000m -jar /usr/gitc/picard.jar \
            CollectRawWgsMetrics \
            INPUT=${input_bam} \
            VALIDATION_STRINGENCY=SILENT \
            REFERENCE_SEQUENCE=${ref_fasta} \
            INCLUDE_BQ_HISTOGRAM=true \
            INTERVALS=${wgs_coverage_interval_list} \
            OUTPUT=${output_prefix} \
            USE_FAST_ALGORITHM=true \
            READ_LENGTH=${read_length}
    }
    runtime {
        docker: docker
        preemptible: preemptible_tries
        memory: "30 GB"
        disks: "local-disk " + disk_size + " HDD"
        maxRetries: max_retries
    }
    output {
       File metrics = "${output_prefix}"
    }
}


task ValidateSamFile {
    File input_bam
    File input_bam_index
    String output_prefix
    File ref_dict
    File ref_fasta
    File ref_fasta_index

    String docker
    Int preemptible_tries
    Int max_retries
    Int disk_size
    Int mem_multiplier = 1
    Int mem_gb = 8 * mem_multiplier

    command {
        java -jar /usr/gitc/picard.jar \
            ValidateSamFile \
            INPUT=${input_bam} \
            OUTPUT=${output_prefix} \
            REFERENCE_SEQUENCE=${ref_fasta} \
            IGNORE="MISSING_TAG_NM" \
            IGNORE="INVALID_VERSION_NUMBER" \
            MODE=VERBOSE \
            IS_BISULFITE_SEQUENCED=false
  }
  runtime {
      docker: docker
      preemptible: preemptible_tries
      memory: "${mem_gb} GB"
      disks: "local-disk " + disk_size + " HDD"
      maxRetries: max_retries
  }
  output {
      File report = "${output_prefix}"
  }
}
