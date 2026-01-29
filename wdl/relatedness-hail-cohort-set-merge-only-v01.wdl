version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/relatedness-hail-v01.wdl" as relatednessHail
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/relatedness-hail-subset-samples-v01.wdl" as relatednessHailSubsetSamples
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/mergeVCFSamples.wdl" as mergeVCFSamples

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow RelatednessCohortSet {
    input {
        Array[File]? somalier_vcf_files
        File? somalier_vcf_file_
        Array[File] ped_uri
        Array[String] cohort_prefixes
        File bed_file
        String merged_filename  # no file extension
        String relatedness_qc_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/scripts/hail_relatedness_check_v0.1.py"
        String plot_relatedness_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/scripts/hail_relatedness_plot_v0.1.py"
        String sex_qc_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/scripts/hail_impute_sex_v0.1.py"
        String sv_base_mini_docker
        String hail_docker
        String bucket_id
        String genome_build
        String x_metric='ibd0'
        String y_metric='kin'
        String kinship_field='kin'  # for sorting in removeDuplicates
        Int chunk_size=100000
        Int samples_per_chunk=0
        Boolean sort_after_merge=false
        Boolean split_multi=true
        RuntimeAttr? runtime_attr_rename_vcf
        RuntimeAttr? runtime_attr_subset_vcfs
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_impute_sex
        RuntimeAttr? runtime_attr_hail_pca
        RuntimeAttr? runtime_attr_check_relatedness
        RuntimeAttr? runtime_attr_plot_relatedness
        RuntimeAttr? runtime_attr_merge_results
    }

    if (!defined(somalier_vcf_file_)) {
        scatter (pair in zip(cohort_prefixes, select_first([somalier_vcf_files]))) {
            call mergeVCFSamples.renameVCFSamplesWithCohort as renameVCFSamplesWithCohort {
                input:
                    vcf_uri=pair.right,
                    cohort_prefix=pair.left,
                    hail_docker=hail_docker,
                    runtime_attr_override=runtime_attr_rename_vcf
            }
        }
        call mergeVCFSamples.mergeVCFSamples as mergeVCFSamples {
            input:
            vcf_files=renameVCFSamplesWithCohort.renamed_vcf_file,
            output_vcf_name=merged_filename+'.vcf.gz',
            sv_base_mini_docker=sv_base_mini_docker,
            format_fields_to_keep=["GT","AD","DP"],
            recalculate_af=true,
            runtime_attr_override=runtime_attr_merge_vcfs
        }
    }

    File merged_vcf_file = select_first([somalier_vcf_file_, mergeVCFSamples.merged_vcf_file])

    call mergePeds {
        input:
            ped_uris=ped_uri,
            cohort_prefixes=cohort_prefixes,
            merged_filename=merged_filename,
            hail_docker=hail_docker,
            input_size=size(ped_uri, 'GB')
    }

    output {
        String cohort_prefix = merged_filename
        File somalier_vcf_file = merged_vcf_file
    }
}

task mergePeds {
     input {
        Array[String] cohort_prefixes
        Array[String] ped_uris
        String hail_docker
        String merged_filename
        Float input_size
        RuntimeAttr? runtime_attr_override
    }

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
        cat <<EOF > merge_peds.py
        import pandas as pd
        import numpy as np
        import sys

        tsvs = pd.read_csv(sys.argv[1], header=None)[0].tolist()
        cohort_prefixes = pd.read_csv(sys.argv[2], header=None)[0].tolist()
        merged_filename = sys.argv[3] + '.ped'

        dfs = []
        tot = len(tsvs)
        for i, (cohort, tsv) in enumerate(zip(cohort_prefixes, tsvs)):
            print(f"Loading pedigree {i+1}/{tot}...")
            df = pd.read_csv(tsv, sep='\t').iloc[:, :6]
            df.columns = ['FamilyID', 'IndividualID', 'FatherID', 'MotherID', 'Sex', 'Affected']
            df[['FamilyID', 'IndividualID', 'FatherID', 'MotherID']] = df[['FamilyID', 'IndividualID', 'FatherID', 'MotherID']].astype(str)
            df[['FamilyID', 'IndividualID']] = f"{cohort}:" + df[['FamilyID', 'IndividualID']]

            df['FatherID'] = df.FatherID.map({samp: f"{cohort}:" + samp for samp in df[df.FatherID!='0'].FatherID.tolist()} | {'0': '0'})
            df['MotherID'] = df.MotherID.map({samp: f"{cohort}:" + samp for samp in df[df.MotherID!='0'].MotherID.tolist()} | {'0': '0'})

            dfs.append(df)
        merged_df = pd.concat(dfs)
        merged_df.to_csv(merged_filename, sep='\t', index=False)
        EOF

        python3 merge_peds.py ~{write_lines(ped_uris)} ~{write_lines(cohort_prefixes)} ~{merged_filename} > stdout
    >>>

    output {
        File merged_ped_file = "~{merged_filename}.ped" 
    }   
}