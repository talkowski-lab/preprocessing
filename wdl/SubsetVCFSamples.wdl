version 1.0

import "helpers.wdl" as helpers
import "mergeVCFs.wdl" as mergeVCFs
import "mergeVCFSamples.wdl" as mergeVCFSamples

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow SubsetVCFSamples {
    input {
        Array[File] vcf_files
        Array[String] samples_array

        String sv_base_mini_docker
    }

    scatter (vcf_file in vcf_files) {
        call helpers.subsetVCFSamples as subsetVCFSamples {
            input:
            vcf_file=vcf_file,
            samples_file=write_lines(samples_array),
            docker=sv_base_mini_docker
        }
    }

    output {
        Array[File] output_vcf_files = subsetVCFSamples.vcf_subset
        Array[File] output_vcf_idxs = subsetVCFSamples.vcf_subset_idx
    }
}