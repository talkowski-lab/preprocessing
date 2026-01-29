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
        Boolean recalculate_af=false
        Boolean fill_missing=false
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
        if [ "~{fill_missing}" = "true" ]; then
            # move merged VCF to temporary file so final output has same filename, regardless of recalculate_af
            mv ~{output_vcf_name} tmp_~{output_vcf_name}
            bcftools +setGT -Oz -o ~{output_vcf_name} tmp_~{output_vcf_name} -- -t . -n 0
        fi
        if [ "~{recalculate_af}" = "true" ]; then
            # move merged VCF to temporary file so final output has same filename, regardless of recalculate_af
            mv ~{output_vcf_name} tmp_~{output_vcf_name}
            bcftools +fill-tags tmp_~{output_vcf_name} -Oz -o ~{output_vcf_name} -- -t AN,AC,AF
        fi
        tabix ~{output_vcf_name}
    >>>

    output {
        File merged_vcf_file = output_vcf_name
        File merged_vcf_idx = output_vcf_name + '.tbi'
    }
}

task mergeVCFSamples {
    input {
        Array[File] vcf_files
        String output_vcf_name
        String sv_base_mini_docker

        # NEW: FORMAT fields to keep; default = keep all
        Array[String]? format_fields_to_keep

        Boolean recalculate_af=false
        Boolean fill_missing=false
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

    # Build FORMAT exclusion expression if requested
    FORMAT_EXCLUDE=""
    if [ "~{defined(format_fields_to_keep)}" = "true" ]; then
        if [ "~{length(select_first([format_fields_to_keep]))}" -gt 0 ]; then
            keep_fields="~{sep=',' format_fields_to_keep}"
            FORMAT_EXCLUDE="-x ^FORMAT/${keep_fields//,/\,FORMAT/}"
        fi
    fi

    processed_vcfs=""
    for vcf in $(cat vcfs_sorted.list); do
        base=$(basename $vcf .vcf.gz)
        processed_vcf=${base}.proc.vcf.gz

        cmd="bcftools annotate"

        # Drop AF if recalculating
        if [ "~{recalculate_af}" = "true" ]; then
            cmd="$cmd -x INFO/AF,FORMAT/AF"
        fi

        # Apply FORMAT field filtering if specified
        if [ -n "$FORMAT_EXCLUDE" ]; then
            cmd="$cmd $FORMAT_EXCLUDE"
        fi

        # If no modifications, symlink original
        if [ "$cmd" = "bcftools annotate" ]; then
            ln -s $vcf $processed_vcf
            tabix -f $processed_vcf
        else
            $cmd -Oz -o $processed_vcf $vcf
            tabix $processed_vcf
        fi

        processed_vcfs="$processed_vcfs $processed_vcf"
    done

    printf "%s\n" $processed_vcfs > vcfs_use.list

    bcftools merge -m none -Oz -o ~{output_vcf_name} --file-list vcfs_use.list

    if [ "~{fill_missing}" = "true" ]; then
        mv ~{output_vcf_name} tmp_~{output_vcf_name}
        bcftools +setGT -Oz -o ~{output_vcf_name} tmp_~{output_vcf_name} -- -t . -n 0
    fi

    if [ "~{recalculate_af}" = "true" ]; then
        mv ~{output_vcf_name} tmp_~{output_vcf_name}
        bcftools +fill-tags tmp_~{output_vcf_name} -Oz -o ~{output_vcf_name} -- -t AN,AC,AF
    fi

    tabix ~{output_vcf_name}
    >>>

    output {
        File merged_vcf_file = output_vcf_name
        File merged_vcf_idx = output_vcf_name + ".tbi"
    }
}

# task mergeVCFSamples {
#     input {
#         Array[File] vcf_files
#         String output_vcf_name
#         String sv_base_mini_docker
#         Boolean keep_gt_ad_dp_only=false
#         Boolean recalculate_af=false
#         Boolean fill_missing=false
#         RuntimeAttr? runtime_attr_override
#     }

#     Float input_size = size(vcf_files, 'GB')
#     Float base_disk_gb = 10.0
#     Float input_disk_scale = 5.0

#     RuntimeAttr runtime_default = object {
#         mem_gb: 4,
#         disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
#         cpu_cores: 1,
#         preemptible_tries: 3,
#         max_retries: 1,
#         boot_disk_gb: 10
#     }

#     RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

#     Float memory = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
#     Int cpu_cores = select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    
#     runtime {
#         memory: "~{memory} GB"
#         disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
#         cpu: cpu_cores
#         preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
#         maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
#         docker: sv_base_mini_docker
#         bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
#     }

#     command <<<
#     set -euo pipefail
#     VCFS="~{write_lines(vcf_files)}"
#     cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list

#     # Preprocess VCFs depending on options
#     processed_vcfs=""
#     for vcf in $(cat vcfs_sorted.list); do
#         base=$(basename $vcf .vcf.gz)
#         processed_vcf=${base}.proc.vcf.gz

#         cmd="bcftools annotate"

#         # If recalculate_af, drop AF fields
#         if [ "~{recalculate_af}" = "true" ]; then
#             cmd="$cmd -x INFO/AF,FORMAT/AF"
#         fi

#         # If keep_gt_ad_dp_only, drop all FORMAT fields except GT
#         if [ "~{keep_gt_ad_dp_only}" = "true" ]; then
#             cmd="$cmd -x ^FORMAT/GT,FORMAT/AD,FORMAT/DP"
#         fi

#         # If no modifications, symlink original
#         if [ "$cmd" = "bcftools annotate" ]; then
#             ln -s $vcf $processed_vcf
#             tabix -f $processed_vcf
#         else
#             $cmd -Oz -o $processed_vcf $vcf
#             tabix $processed_vcf
#         fi

#         processed_vcfs="$processed_vcfs $processed_vcf"
#     done

#     # Write processed VCFs to new list
#     printf "%s\n" $processed_vcfs > vcfs_use.list

#     # Merge processed VCFs
#     bcftools merge -m none -Oz -o ~{output_vcf_name} --file-list vcfs_use.list

#     if [ "~{fill_missing}" = "true" ]; then
#         mv ~{output_vcf_name} tmp_~{output_vcf_name}
#         bcftools +setGT -Oz -o ~{output_vcf_name} tmp_~{output_vcf_name} -- -t . -n 0
#     fi

#     if [ "~{recalculate_af}" = "true" ]; then
#         mv ~{output_vcf_name} tmp_~{output_vcf_name}
#         bcftools +fill-tags tmp_~{output_vcf_name} -Oz -o ~{output_vcf_name} -- -t AN,AC,AF
#     fi

#     tabix ~{output_vcf_name}
#     >>>

#     output {
#         File merged_vcf_file = output_vcf_name
#         File merged_vcf_idx = output_vcf_name + '.tbi'
#     }
# }

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
            bcftools +fill-tags tmp_~{output_vcf_name} -Oz -o ~{output_vcf_name} -- -t AF
        fi
        tabix ~{output_vcf_name}
    >>>

    output {
        File merged_vcf_file = output_vcf_name
        File merged_vcf_idx = output_vcf_name + '.tbi'
    }
}

task renameVCFSamplesWithCohort {
    input {
        File vcf_uri
        String cohort_prefix
        String hail_docker
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_uri, 'GB')
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
    cat <<EOF > merge_vcfs.py
    import pandas as pd
    import numpy as np
    import hail as hl
    import sys
    import os

    vcf_uri = sys.argv[1]
    cohort_prefix = sys.argv[2]
    cores = sys.argv[3]  # string
    mem = int(np.floor(float(sys.argv[4])))

    hl.init(min_block_size=128, spark_conf={"spark.executor.cores": cores, 
                        "spark.executor.memory": f"{int(np.floor(mem*0.4))}g",
                        "spark.driver.cores": cores,
                        "spark.driver.memory": f"{int(np.floor(mem*0.4))}g"
                        }, tmp_dir="tmp", local_tmpdir="tmp")

    mt = hl.import_vcf(vcf_uri, reference_genome='GRCh38', force_bgz=True, call_fields=[], array_elements_required=False)
    mt = mt.key_cols_by()
    mt = mt.annotate_cols(s=hl.str(f"{cohort_prefix}:")+mt.s).key_cols_by('s')
    try:
        # for haploid (e.g. chrY)
        mt = mt.annotate_entries(
            GT = hl.if_else(
                    mt.GT.ploidy == 1, 
                    hl.call(mt.GT[0], mt.GT[0]),
                    mt.GT)
        )
    except:
        pass

    header = hl.get_vcf_metadata(vcf_uri) 
    hl.export_vcf(mt, f"{os.path.basename(vcf_uri).split('.vcf')[0]}_new_sample_IDs.vcf.bgz", metadata=header)
    EOF

    python3 merge_vcfs.py ~{vcf_uri} ~{cohort_prefix} ~{cpu_cores} ~{memory}
    >>>

    String file_ext = if sub(basename(vcf_uri), '.vcf.gz', '')!=basename(vcf_uri) then '.vcf.gz' else '.vcf.bgz'

    output {
        File renamed_vcf_file = "~{basename(vcf_uri, file_ext)}_new_sample_IDs.vcf.bgz"
    }
}
