version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/mergeVCFSamples.wdl" as mergeVCFSamples
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/helpers.wdl" as helpers
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/wdl/ancestry-inference-hail-v01.wdl" as AncestryInference

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}

workflow AncestryInferenceCohortSet {
    input {
        # Array[Array[File]] vep_vcf_files
        # Array[String] cohort_prefixes

        Array[File]? ancestry_vcf_files
        File? ancestry_vcf_file_
        File gnomad_vcf_uri
        File gnomad_rf_onnx
        File pop_labels_tsv
        String gnomad_loading_ht
        String infer_ancestry_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/eren_dev/scripts/hail_infer_ancestry_v0.1.py"

        Array[String] cohort_prefixes
        String cohort_prefix
        String hail_docker
        String sv_base_mini_docker

        Int num_pcs=16  # for gnomADv3
        Float min_prob=0.75

        Boolean infer_ancestry=true
        Boolean use_gnomad_rf=false
        Boolean filter_entries_before=true

        String genome_build='GRCh38'

        RuntimeAttr? runtime_attr_merge_vcfs
        RuntimeAttr? runtime_attr_infer_ancestry
    }

    if ((!defined(ancestry_vcf_file_)) || (ancestry_vcf_file_ == '')) {        
        scatter (pair in zip(cohort_prefixes, select_first([ancestry_vcf_files]))) {
            call mergeVCFSamples.renameVCFSamplesWithCohort as renameVCFSamples {
                input:
                    vcf_uri=pair.right,
                    cohort_prefix=pair.left,
                    hail_docker=hail_docker
            }
        }
        call mergeVCFSamples.mergeVCFSamples as mergeVCFSamples {
            input:
            vcf_files=renameVCFSamples.renamed_vcf_file,
            output_vcf_name=cohort_prefix+'.ancestry_sites.vcf.gz',
            sv_base_mini_docker=sv_base_mini_docker,
            format_fields_to_keep=["GT","GQ","DP","AD"],
            recalculate_af=true,
            runtime_attr_override=runtime_attr_merge_vcfs
        }
    }

    File vcf_uri = select_first([ancestry_vcf_file_, mergeVCFSamples.merged_vcf_file])

    if (infer_ancestry) {
        call AncestryInference.inferAncestry as inferAncestry {
            input:
                vcf_uri=vcf_uri,
                gnomad_vcf_uri=gnomad_vcf_uri,
                gnomad_rf_onnx=gnomad_rf_onnx,
                pop_labels_tsv=pop_labels_tsv,
                cohort_prefix=cohort_prefix,
                gnomad_loading_ht=gnomad_loading_ht,
                infer_ancestry_script=infer_ancestry_script,
                hail_docker=hail_docker,
                num_pcs=num_pcs,
                min_prob=min_prob,
                use_gnomad_rf=use_gnomad_rf,
                filter_entries_before=filter_entries_before,
                genome_build=genome_build,
                runtime_attr_override=runtime_attr_infer_ancestry
        }
    }
    output {
        File ancestry_vcf_file = vcf_uri
        File ancestry_tsv = select_first([inferAncestry.ancestry_tsv])
        File ancestry_plot = select_first([inferAncestry.ancestry_plot])
    }
}