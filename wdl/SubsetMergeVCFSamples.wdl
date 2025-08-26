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

    scatter (vcf_file in flatten(cohort_vcf_files)) {
        call helpers.subsetVCFSamples as subsetCohortVCFSamples {
            input:
            vcf_file=vcf_file,
            samples_file=write_lines(samples_array),
            docker=sv_base_mini_docker
        }
    }

    call mergeVCFSamples.mergeVCFs as mergeVCFSamples {
        input:
        vcf_files=subsetCohortVCFSamples.vcf_subset,
        output_vcf_name=merged_prefix + '.merged.vcf.gz',
        sv_base_mini_docker=sv_base_mini_docker
    }

    output {
        File merged_vcf_file = mergeVCFSamples.merged_vcf_file
        File merged_vcf_idx = mergeVCFSamples.merged_vcf_idx
    }
}