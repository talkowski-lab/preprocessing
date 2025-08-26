version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/helpers.wdl" as helpers
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/mergeVCFs.wdl" as mergeVCFs
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/mergeVCFSamples.wdl" as mergeVCFSamples

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow MergeVCFSamples {
    input {
        Array[File] cohort_vcf_files
        Array[String] samples_array

        Array[String] cohort_prefixes
        String merged_prefix
        String sv_base_mini_docker
    }

    scatter (vcf_file in cohort_vcf_files) {
        call helpers.subsetVCFSamples as subsetCohortVCFSamples {
            input:
            vcf_file=vcf_file,
            samples_file=write_lines(samples_array),
            docker=sv_base_mini_docker
        }
    }

    call mergeCommonVCFs as mergeCommonVCFs {
        input:
        vcf_files=subsetCohortVCFSamples.vcf_subset,
        output_vcf_name=merged_prefix + '.merged.vcf.gz',
        sv_base_mini_docker=sv_base_mini_docker
    }

    output {
        File merged_vcf_file = mergeCommonVCFs.merged_vcf_file
        File merged_vcf_idx = mergeCommonVCFs.merged_vcf_idx
    }
}

task mergeCommonVCFs {
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

        # Index inputs
        for vcf in $(cat vcfs_sorted.list); do
            tabix $vcf
        done

        # Get number of VCFs
        N=$(wc -l < vcfs_sorted.list)

        # Step 1: Extract only sites common across all VCFs
        bcftools isec -n =$N -w1 -Oz -o common_sites.vcf.gz $(cat vcfs_sorted.list)
        tabix common_sites.vcf.gz

        # Step 2: Merge VCFs, restricted to common sites
        bcftools merge -m none -Oz -o ~{output_vcf_name} --file-list vcfs_sorted.list --no-update -R common_sites.vcf.gz

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
