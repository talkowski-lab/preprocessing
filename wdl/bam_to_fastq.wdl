version 1.0

workflow BamToFastq {
  input {
    File input_reads
    String output_prefix

    File? ref_fasta
    File? ref_fasta_fai
    Int bam_to_fastq_disk_gb
  }

  call SamtoolsBamToFastq {
    input:
      input_reads   = input_reads,
      output_prefix = output_prefix,
      ref_fasta     = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      bam_to_fastq_disk_gb = bam_to_fastq_disk_gb
  }

  output {
    File fastq_1 = SamtoolsBamToFastq.fastq_1
    File fastq_2 = SamtoolsBamToFastq.fastq_2
  }
}

task SamtoolsBamToFastq {
  input {
    File input_reads
    String output_prefix
    File? ref_fasta
    File? ref_fasta_fai
    Int bam_to_fastq_disk_gb
  }

  Boolean is_cram    = sub(input_reads, ".*\\.", "") == "cram"
  Boolean has_ref    = defined(ref_fasta)

  command <<<
    set -euo pipefail

    # 1. Determine if a reference is required (for CRAMs)
    REF_ARGS=""
    if ~{true="true" false="false" is_cram} ; then
      if ! ~{true="true" false="false" has_ref} ; then
        echo "Error: Reference FASTA is required for CRAM files." >&2
        exit 1
      fi
      REF_ARGS="--reference ~{ref_fasta}"
    fi

    # 2. Collate alignments by queryname
    samtools collate -@ 4 $REF_ARGS -o collated.bam -T tmp_collate "~{input_reads}"

    # 3. Convert to paired compressed FASTQs
    samtools fastq -@ 4 \
      -1 "~{output_prefix}_1.fastq.gz" \
      -2 "~{output_prefix}_2.fastq.gz" \
      -0 /dev/null -s /dev/null \
      collated.bam
  >>>

  runtime {
    docker:      "staphb/samtools:1.17"
    cpu:         4
    memory:      "16 GB"
    disks:       "local-disk " + bam_to_fastq_disk_gb + " SSD"
    preemptible: 1
  }

  output {
    File fastq_1 = "~{output_prefix}_1.fastq.gz"
    File fastq_2 = "~{output_prefix}_2.fastq.gz"
  }
}
