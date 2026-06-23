version 1.0

workflow AddDepthFields {
    input {
        Array[File] sharded_vcfs
        Array[File] sharded_vcf_indices
        String output_prefix
        File? reference_fasta
        File? reference_index
    }

    scatter (i in range(length(sharded_vcfs))) {
        call AddDepthToVCF {
            input:
                vcf = sharded_vcfs[i],
                vcf_index = sharded_vcf_indices[i],
                shard_name = basename(sharded_vcfs[i], ".vcf.gz"),
                reference_fasta = reference_fasta,
                reference_index = reference_index
        }
    }

    output {
        Array[File] vcfs_with_depth = AddDepthToVCF.output_vcf
        Array[File] vcf_indices = AddDepthToVCF.output_vcf_index
    }
}

task AddDepthToVCF {
    input {
        File vcf
        File vcf_index
        String shard_name
        File? reference_fasta
        File? reference_index
        
        Int disk_size = ceil(size(vcf, "GB") * 3) + 20
        Int mem_gb = 16
        Int cpu = 2
    }

    command <<<
        set -euo pipefail

        # Install tabix
        apt-get update && apt-get install -y tabix

        # Create Python script to run Hail
        cat << 'EOF' > add_depth.py
import hail as hl
import sys

input_vcf = sys.argv[1]
output_vcf = sys.argv[2]

hl.init(log="/dev/null")

# Import VCF with GRCh38 reference genome and allow missing array elements
mt = hl.import_vcf(input_vcf, force_bgz=True, reference_genome='GRCh38', array_elements_required=False)

# Add FORMAT/DP = sum(AD)
mt = mt.annotate_entries(
    DP = hl.sum(mt.AD)
)

# Add INFO/DP = sum of sample DP values (add to info struct)
mt = mt.annotate_rows(
    info = mt.info.annotate(DP = hl.agg.sum(mt.DP))
)

# Export final VCF
header = hl.get_vcf_metadata(input_vcf)
hl.export_vcf(mt, output_vcf, metadata=header)
EOF

        # Run hail script
        python3 add_depth.py "~{vcf}" "~{shard_name}.with_depth.vcf.bgz"

        # Rename and index the VCF
        mv "~{shard_name}.with_depth.vcf.bgz" "~{shard_name}.with_depth.vcf.gz"
        tabix -p vcf "~{shard_name}.with_depth.vcf.gz"
    >>>

    output {
        File output_vcf = "~{shard_name}.with_depth.vcf.gz"
        File output_vcf_index = "~{shard_name}.with_depth.vcf.gz.tbi"
    }

    runtime {
        docker: "us.gcr.io/talkowski-sv-gnomad/shineren:hail-0.2.134"
        memory: "~{mem_gb} GB"
        disks: "local-disk ~{disk_size} HDD"
        cpu: cpu
        preemptible: 3
    }
}
