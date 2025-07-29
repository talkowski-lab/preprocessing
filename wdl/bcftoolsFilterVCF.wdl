version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow BcftoolsFilter {
    input {
        File vcf_file
        String filter_command
        String output_filename
        String docker
    }

    call bcftoolsFilter {
        input:
        vcf_file=vcf_file,
        filter_command=filter_command,
        output_filename=output_filename,
        docker=docker
    }

    output {
        File output_vcf = bcftoolsFilter.output_vcf
        File output_vcf_idx = bcftoolsFilter.output_vcf_idx
    }
}

task bcftoolsFilter {
    input {
        File vcf_file
        String filter_command
        String output_filename
        String docker = "us.gcr.io/prenatal-wgs-sv/cwhelan/cffdna-pyro:bcftools221_391ecb"
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_file, "GB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 10.0
    RuntimeAttr runtime_default = object {
        mem_gb: 8,
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
        docker: docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command {
        set -euo pipefail

        ~{filter_command} ~{vcf_file} -Oz -o ~{output_filename}
        tabix ~{output_filename}
    }

    output {
        File output_vcf = output_filename
        File output_vcf_idx = output_filename + '.tbi'
    }
}