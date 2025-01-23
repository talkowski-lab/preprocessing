version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow MergeVCFs {
    input {
        # give either text file with VCF paths or Array[File] of VCF paths
        File? vcf_list_file
        Array[File]? vcf_files
        # if missing tags in header gives an error with bcftools merge, mainly for DDD (probably safe to ignore this parameter for other cohorts)
        File? header_file
        # hacky renaming, just for DDD (safe to ignore this parameter for other cohorts)
        Boolean rename_samples=false
        # recalculate AF from AC and AN after merge
        Boolean recalculate_af=false
        String sample_set_id
        String sv_base_mini_docker
    }

    if (rename_samples) {
        call renameVCFSamples {
            input:
            vcf_files=select_first([vcf_files, read_lines(select_first([vcf_list_file]))]),
            sv_base_mini_docker=sv_base_mini_docker
        }
    }

    if (defined(header_file)) {
        call mergeVCFsReheader {
            input:
            vcf_files=select_first([renameVCFSamples.renamed_vcf_files, vcf_files, read_lines(select_first([vcf_list_file]))]),
            header_file=select_first([header_file]),
            output_vcf_name=sample_set_id + '.merged.vcf.gz',
            recalculate_af=recalculate_af,
            sv_base_mini_docker=sv_base_mini_docker
        }
    }

    if (!defined(header_file)) {
        call mergeVCFs {
            input:
            vcf_files=select_first([renameVCFSamples.renamed_vcf_files, vcf_files, read_lines(select_first([vcf_list_file]))]),
            output_vcf_name=sample_set_id + '.merged.vcf.gz',
            sv_base_mini_docker=sv_base_mini_docker,
            recalculate_af=recalculate_af
        }
    }

    output {
        File merged_vcf_file = select_first([mergeVCFs.merged_vcf_file, mergeVCFsReheader.merged_vcf_file])
        File merged_vcf_idx = select_first([mergeVCFs.merged_vcf_idx, mergeVCFsReheader.merged_vcf_idx])
    }
}

task renameVCFSamples {
    input {
        Array[File] vcf_files
        String sv_base_mini_docker
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        for vcf in $(cat vcfs_sorted.list);
        do
            bcftools query -l $vcf > samples.txt;
            trio_id=$(basename $vcf | cut -f 1 -d '.');
            sed "s/^/$trio_id\-/" samples.txt > new_samples.txt;
            paste samples.txt new_samples.txt > sample_map.txt;
            bcftools reheader -s sample_map.txt -o $(basename $vcf).renamed.vcf.gz $vcf;
        done
    >>>

    output {
        Array[File] renamed_vcf_files = glob('*.renamed.vcf.gz')
    }
}

task mergeVCFs {
    input {
        Array[File] vcf_files
        String output_vcf_name
        String sv_base_mini_docker
        Boolean recalculate_af
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        for vcf in $(cat vcfs_sorted.list);
        do
            tabix $vcf;
        done
        bcftools merge -m none -Oz -o ~{output_vcf_name} --file-list vcfs_sorted.list
        if [ "~{recalculate_af}" = "true" ]; then
            # move merged VCF to temporary file so final output has same filename, regardless of recalculate_af
            mv ~{output_vcf_name} tmp_~{output_vcf_name}
            bcftools +fill-tags tmp_~{output_vcf_name} -Oz -o ~{output_vcf_name} -- -t INFO/AF
        fi
        tabix ~{output_vcf_name}
    >>>

    output {
        File merged_vcf_file = output_vcf_name
        File merged_vcf_idx = output_vcf_name + '.tbi'
    }
}

task mergeVCFsReheader {
    input {
        Array[File] vcf_files
        File header_file
        String output_vcf_name
        String sv_base_mini_docker
        Boolean recalculate_af
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        VCFS="~{write_lines(vcf_files)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        for vcf in $(cat vcfs_sorted.list);
        do
            bcftools annotate -h ~{header_file} -Oz -o "$vcf".reheader.vcf.gz $vcf
            tabix "$vcf".reheader.vcf.gz
            echo "$vcf".reheader.vcf.gz >> vcfs_sorted_reheader.list;
        done
        bcftools merge -m id -Oz -o ~{output_vcf_name} --file-list vcfs_sorted_reheader.list
        if [ "~{recalculate_af}" = "true" ]; then
            # move merged VCF to temporary file so final output has same filename, regardless of recalculate_af
            mv ~{output_vcf_name} tmp_~{output_vcf_name}
            bcftools +fill-tags tmp_~{output_vcf_name} -Oz -o ~{output_vcf_name} -- -t INFO/AF
        fi
        tabix ~{output_vcf_name}
    >>>

    output {
        File merged_vcf_file = output_vcf_name
        File merged_vcf_idx = output_vcf_name + '.tbi'
    }
}