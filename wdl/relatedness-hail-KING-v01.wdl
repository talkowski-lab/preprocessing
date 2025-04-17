version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/main/wdl/mergeVCFs.wdl" as mergeVCFs
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/main/wdl/helpers.wdl" as helpers
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/relatedness-hail-v01.wdl" as relatednessHail

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow Relatedness {
    input {
        Array[File]? vep_vcf_files
        File? somalier_vcf_file_
        File ped_uri
        File bed_file
        String cohort_prefix
        String relatedness_qc_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/scripts/hail_relatedness_check_KING_v0.1.py"
        String plot_relatedness_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/scripts/hail_relatedness_plot_v0.1.py"
        String sv_base_mini_docker
        String hail_docker
        String bucket_id
        String genome_build
        String x_metric='ibs0'
        String y_metric='phi'
        Int chunk_size=0
        Boolean sort_after_merge=false
        Boolean split_multi=true
        RuntimeAttr? runtime_attr_subset_vcfs
        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_impute_sex
        RuntimeAttr? runtime_attr_check_relatedness
        RuntimeAttr? runtime_attr_plot_relatedness
    }

    if (!defined(somalier_vcf_file_)) {
        
        scatter (vcf_uri in select_first([vep_vcf_files])) {
            String filename = basename(vcf_uri)
            String prefix = if (sub(filename, ".gz", "")!=filename) then basename(filename, ".vcf.gz") else basename(filename, ".vcf.bgz")
            call helpers.subsetVCFs as subsetVCFs {
                input:
                    bed_file=bed_file,
                    vcf_uri=vcf_uri,
                    vcf_idx=vcf_uri+'.tbi',
                    output_name=prefix + '.somalier.subset.vcf.gz',
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_subset_vcfs
            }
        }

        call mergeVCFs.mergeVCFs as mergeVCFs {
        input:
            vcf_files=subsetVCFs.subset_vcf,
            sv_base_mini_docker=sv_base_mini_docker,
            cohort_prefix=cohort_prefix,
            sort_after_merge=sort_after_merge,
            runtime_attr_override=runtime_attr_merge_vcfs
        }
    }
    File merged_vcf_file = select_first([somalier_vcf_file_, mergeVCFs.merged_vcf_file])

    call relatednessHail.checkRelatedness {
        input:
        vcf_uri=merged_vcf_file,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        relatedness_qc_script=relatedness_qc_script,
        hail_docker=hail_docker,
        bucket_id=bucket_id,
        genome_build=genome_build,
        runtime_attr_override=runtime_attr_check_relatedness
    }

    call relatednessHail.plotRelatedness {
        input:
        kinship_tsv=checkRelatedness.kinship_tsv,
        ped_uri=ped_uri,
        cohort_prefix=cohort_prefix,
        plot_relatedness_script=plot_relatedness_script,
        hail_docker=hail_docker,
        x_metric=x_metric,
        y_metric=y_metric,
        chunk_size=chunk_size,
        runtime_attr_override=runtime_attr_plot_relatedness
    }

    output {
        File somalier_vcf_file = merged_vcf_file
        File relatedness_qc = checkRelatedness.relatedness_qc
        File relatedness_plot = plotRelatedness.relatedness_plot
    }
}