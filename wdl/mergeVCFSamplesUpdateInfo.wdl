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
        String merge_vcf_samples_update_info_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/scripts/hail_merge_vcf_samples_update_info.py"
        String hail_docker
        String prefix

        Int batch_size = 10
        String? output_vcf_filename_override  # End with ".bgz" for bgzipped output
        Array[String] min_info_fields = ["FractionInformativeReads", "LOD", "MQ", "MQRankSum", "BaseQRankSum", "ReadPosRankSum", "QD", "VQSLOD", "QUALapprox", "ExcessHet"]
        Array[String] max_info_fields = ["FS", "SOR"]
        Array[String] sum_info_fields = ["DP", "MQ_DP", "VarDP"]
    }

    String output_vcf_file = select_first([output_vcf_filename_override, prefix + '.merged.vcf.bgz'])

    call mergeVCFSamplesUpdateInfo {
        input:
            merge_vcf_samples_update_info_script = merge_vcf_samples_update_info_script,
            vcf_files = vcf_files,
            batch_size = batch_size,
            output_vcf_file = output_vcf_file,
            min_info_fields = min_info_fields,
            max_info_fields = max_info_fields,
            sum_info_fields = sum_info_fields,
            hail_docker = hail_docker
    }

    output {
        File merged_vcf_file = mergeVCFSamplesUpdateInfo.merged_vcf_file
    }
}

task mergeVCFSamplesUpdateInfo {
    input {
        String merge_vcf_samples_update_info_script
        Array[File] vcf_files
        Int batch_size
        String output_vcf_file
        String hail_docker
        Array[String] min_info_fields
        Array[String] max_info_fields
        Array[String] sum_info_fields
        RuntimeAttr? runtime_attr_override
    }   

    Float input_size = size(vcf_files, 'GB')
    Float base_disk_gb = 10.0
    Float input_disk_scale = 5.0

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

        # Run the Python script with the provided arguments
        python3 merge_vcf_files.py \
            --vcf-files ~{sep=" " vcf_files} \
            --batch-size ~{batch_size} \
            --output-vcf-file ~{output_vcf_file} \
            --mem ~{memory} \
            --min-info-fields ~{sep=" " min_info_fields} \
            --max-info-fields ~{sep=" " max_info_fields} \
            --sum-info-fields ~{sep=" " sum_info_fields} 
    >>>

    output {
        File merged_vcf_file = output_vcf_file
    }
}
