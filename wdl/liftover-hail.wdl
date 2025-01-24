version 1.0

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

# assumes lifting over from hg19/GRCh37 to hg38/GRCh38 (for now)
workflow liftOver {
    input {
        File vcf_file
        String hail_docker
        String genome_build='GRCh37'
    }

    call liftOverVCF {
        input:
        vcf_file=vcf_file,
        hail_docker=hail_docker,
        genome_build=genome_build
    }

    output {
        File output_vcf_file = liftOverVCF.output_vcf_file
        File output_vcf_idx = liftOverVCF.output_vcf_idx
    }
}

task liftOverVCF {
    input {
        File vcf_file
        String hail_docker
        String genome_build
        RuntimeAttr? runtime_attr_override
    }
    Float input_size = size(vcf_file, "GB")
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

    String file_ext = if sub(basename(vcf_file), '.vcf.gz', '')!=basename(vcf_file) then '.vcf.gz' else '.vcf.bgz'
    String output_filename = basename(vcf_file, file_ext) + '.lifted_over_GRCh38.vcf.bgz'

    command <<<
    set -euo pipefail
    cat <<EOF > liftover.py
    import datetime
    import pandas as pd
    import hail as hl
    import numpy as np
    import sys
    import os
    
    vcf_file = sys.argv[1]
    cores = sys.argv[2]
    mem = int(np.floor(float(sys.argv[3])))
    build = sys.argv[4]
    output_filename = sys.argv[5]

    hl.init(min_block_size=128, 
            local=f"local[*]", 
            spark_conf={
                        "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                        "spark.speculation": 'true'
                        }, 
            tmp_dir="tmp", local_tmpdir="tmp",
                        )
    header = hl.get_vcf_metadata(vcf_file)
    mt = hl.import_vcf(vcf_file, force_bgz=vcf_file.split('.')[-1] in ['gz', 'bgz'], 
                reference_genome=build, array_elements_required=False)
    rg37 = hl.get_reference('GRCh37')  
    rg38 = hl.get_reference('GRCh38')  
    rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)  

    mt = mt.key_rows_by().annotate_rows(locus=hl.liftover(mt.locus, 'GRCh38'))  
    n_rows_before = mt.count_rows()
    mt = mt.filter_rows(hl.is_defined(mt.locus))  
    n_rows_after = mt.count_rows()
    print(f"{n_rows_after}/{n_rows_before} remaining after liftover.")
    mt = mt.key_rows_by('locus','alleles')  
    hl.export_vcf(mt, output_filename, metadata=header, tabix=True)
    EOF 
    python3 liftover.py ~{vcf_file} ~{cpu_cores} ~{memory} ~{genome_build} ~{output_filename}
    >>>

    output {
        File output_vcf_file = output_filename
        File output_vcf_idx = output_filename + '.tbi'
    }
}