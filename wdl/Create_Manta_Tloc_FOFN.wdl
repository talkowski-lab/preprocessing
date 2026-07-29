version 1.0

workflow CreateMantaTlocFOFN {
  input {
    Array[String] manta_tloc_vcfs
    String output_prefix
  }

  call WriteFOFN {
    input:
      vcfs = manta_tloc_vcfs,
      output_prefix = output_prefix
  }

  output {
    File manta_tloc_fofn = WriteFOFN.fofn
  }
}

task WriteFOFN {
  input {
    Array[String] vcfs
    String output_prefix
  }

  command <<<
    set -euo pipefail
    cat ~{write_lines(vcfs)} | sed 's|/mnt/disks/cromwell_root/|gs://|g' > ~{output_prefix}.manta_tloc_vcfs.fofn
  >>>

  output {
    File fofn = "~{output_prefix}.manta_tloc_vcfs.fofn"
  }

  runtime {
    cpu: 1
    memory: "1 GiB"
    disks: "local-disk 10 HDD"
    docker: "ubuntu:latest"
    preemptible: 3
  }
}
