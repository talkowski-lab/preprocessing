version 1.0

import "helpers.wdl" as helpers

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow mergeResultsPython {
     input {
        File tsv_filenames
        String hail_docker
        String merged_filename
        RuntimeAttr? runtime_attr_override
    }

    call helpers.mergeResultsPython as mergeResults {
        input:
        tsvs=read_lines(tsv_filenames),
        hail_docker=hail_docker,
        merged_filename=merged_filename,
        runtime_attr_override=runtime_attr_override
    }

    output {
        File merged_tsv = mergeResults.merged_tsv
    }
}
