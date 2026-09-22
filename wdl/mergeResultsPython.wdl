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
        Array[String] tsvs
        String hail_docker
        String merged_filename
        RuntimeAttr? runtime_attr_override
    }

    call helpers.mergeResultsPython as mergeResults {
        input:
        tsvs=tsvs,
        hail_docker=hail_docker,
        merged_filename=merged_filename,
        input_size=size(tsvs, 'GB'),
        runtime_attr_override=runtime_attr_override
    }

    output {
        File merged_tsv = mergeResults.merged_tsv
    }
}
