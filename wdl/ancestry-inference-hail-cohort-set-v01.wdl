version 1.0

import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/main/wdl/mergeVCFs.wdl" as mergeVCFs
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/main/wdl/helpers.wdl" as helpers
import "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/main/wdl/ancestry-inference-hail-v01.wdl" as AncestryInference

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
        String infer_ancestry_script = "https://raw.githubusercontent.com/talkowski-lab/preprocessing/refs/heads/main/scripts/hail_infer_ancestry_v0.1.py"

        String cohort_set_id
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
    if (!defined(somalier_vcf_file_)) {
        scatter (pair in zip(cohort_prefixes, select_first([somalier_vcf_files]))) {
            call renameVCFSamples {
                input:
                    vcf_uri=pair.right,
                    cohort_prefix=pair.left,
                    hail_docker=hail_docker,
                    runtime_attr_override=runtime_attr_rename_vcf
            }
        }
        call mergeVCFs.mergeVCFs as mergeVCFs {
            input:
            vcf_files=renameVCFSamples.renamed_vcf_file,
            output_vcf_name=merged_filename+'.vcf.gz',
            sv_base_mini_docker=sv_base_mini_docker,
            recalculate_af=false,
            runtime_attr_override=runtime_attr_merge_vcfs
        }
    }

        call mergeVCFSamples {
            input:
            vcf_files=select_first([ancestry_vcf_files]),
            hail_docker=hail_docker,
            prefix=cohort_set_id,
            genome_build=genome_build,
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
                cohort_prefix=cohort_set_id,
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