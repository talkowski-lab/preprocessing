version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow FillMissingGT {
    input {
        File vcf_file
        String sv_base_mini_docker
    }

    call fillMissingGT_indexVCF {
        input:
        vcf_file=vcf_file,
        sv_base_mini_docker=sv_base_mini_docker
    }

    output {
        File output_vcf_file = fillMissingGT_indexVCF.output_vcf_file
        File output_vcf_idx = fillMissingGT_indexVCF.output_vcf_idx
    }
}

task fillMissingGT_indexVCF {
    input {
        File vcf_file
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_override
    }
    
    Float input_size = size([vcf_file], "GB") 
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 4,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
    Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }    

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext)+'.fill.missing.GT.vcf.gz'

    command <<<
    set -euo pipefail
    bcftools filter -S 0 -e 'GT=="./."' -Oz -o ~{output_filename} ~{vcf_file}
    bcftools index -t ~{output_filename}
    >>>

    output {
        File output_vcf_file = output_filename
        File output_vcf_idx = output_filename + '.tbi'
    }
}