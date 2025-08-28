version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow gCNV_BEDtoVCF {
    input {
        File bed_uri

        String gcnv_bed_to_vcf_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/scripts/hail_gcnv_bed_to_vcf.py"
        String hail_docker
        String sv_pipeline_docker

        String genome_build = 'GRCh38'

        Array[String] row_key = ['rsid']
        Array[String] col_key = ['sample_fix']
        Array[String] skip_fields = ['chr','start']  # removed
        Array[String] col_fields = ['sample']
        Array[String] entry_fields = ['name', 'GT', 'CN', 'NP', 'QA', 'QS', 'QSE', 'QSS', 'ploidy', 'strand',
            'ID', 'defragmented', 'sf', 'sc', 'PASS_SAMPLE',
            'PASS_FREQ', 'PASS_QS', 'HIGH_QUALITY']
        Array[String] priority_row_fields = ['END','SVTYPE','SVLEN']

        RuntimeAttr? runtime_attr_hail_bed_to_vcf
        RuntimeAttr? runtime_attr_index_vcf
    }

    call BEDtoVCF {
        input:
            bed_uri=bed_uri,
            hail_docker=hail_docker,
            gcnv_bed_to_vcf_script=gcnv_bed_to_vcf_script,
            row_key=row_key,
            col_key=col_key,
            skip_fields=skip_fields,
            col_fields=col_fields,
            entry_fields=entry_fields,
            priority_row_fields=priority_row_fields,
            genome_build=genome_build,
            runtime_attr_override=runtime_attr_hail_bed_to_vcf
    }

    call fillMissingGT_indexVCF {
        input:
            vcf_file=BEDtoVCF.output_vcf,
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_index_vcf
    }

    output {
        File output_vcf_file = fillMissingGT_indexVCF.output_vcf_file
        File output_vcf_idx = fillMissingGT_indexVCF.output_vcf_idx
    }
}

task BEDtoVCF {
    input {
        File bed_uri

        String genome_build
        String hail_docker
        String gcnv_bed_to_vcf_script

        Array[String] row_key
        Array[String] col_key
        Array[String] skip_fields
        Array[String] col_fields
        Array[String] entry_fields
        Array[String] priority_row_fields

        RuntimeAttr? runtime_attr_override
    }
    
    Float input_size = size([bed_uri], "GB") 
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
        docker: hail_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
    set -euo pipefail
    curl ~{gcnv_bed_to_vcf_script} > convert_bed_to_vcf.py
    python3 convert_bed_to_vcf.py \
        --bed-uri ~{bed_uri} \
        ~{if length(row_key) > 0 then "--row-key " + sep(",", row_key) else ""} \
        ~{if length(col_key) > 0 then "--col-key " + sep(",", col_key) else ""} \
        ~{if length(skip_fields) > 0 then "--skip-fields " + sep(",", skip_fields) else ""} \
        ~{if length(col_fields) > 0 then "--col-fields " + sep(",", col_fields) else ""} \
        ~{if length(entry_fields) > 0 then "--entry-fields " + sep(",", entry_fields) else ""} \
        ~{if length(priority_row_fields) > 0 then "--priority-row-fields " + sep(",", priority_row_fields) else ""} \
        --genome-build ~{genome_build}
    >>>

    String file_ext = if sub(basename(bed_uri), '.bed.gz', '')!=basename(bed_uri) then '.bed.gz' else '.bed.bgz'

    output {
        File output_vcf = basename(bed_uri, file_ext) + '.vcf.bgz'
    }
}

task fillMissingGT_indexVCF {
    input {
        File vcf_file
        String sv_pipeline_docker

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
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }    

    String gz_vcf_file = basename(vcf_file, '.vcf.bgz') + '.vcf.gz'
    String vcf_filled_missing_GT = basename(gz_vcf_file, '.vcf.gz')+'.fill.missing.GT.vcf.gz'

    command <<<
    set -euo pipefail
    mv ~{vcf_file} ~{gz_vcf_file}
    bcftools filter -S 0 -e 'GT=="./."' -Oz -o ~{vcf_filled_missing_GT} ~{gz_vcf_file}
    bcftools index -t ~{vcf_filled_missing_GT}
    >>>

    output {
        File output_vcf_file = vcf_filled_missing_GT
        File output_vcf_idx = vcf_filled_missing_GT + '.tbi'
    }
}