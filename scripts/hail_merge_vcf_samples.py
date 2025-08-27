#!/usr/bin/env python3

import argparse
import numpy as np
import hail as hl

def merge_vcfs(vcf_files):
    # Load first VCF
    mt = hl.import_vcf(vcf_files[0], force_bgz=vcf_files[0].split('.')[-1]=='gz')
    
    # Collect common rows across all VCFs
    common_rows = mt.rows()
    for vcf in vcf_files[1:]:
        other_mt = hl.import_vcf(vcf, force_bgz=vcf.split('.')[-1]=='gz')
        common_rows = common_rows.semi_join(other_mt.rows())

    # Filter the first VCF to common rows
    mt = mt.semi_join_rows(common_rows)

    # Merge the rest
    for vcf in vcf_files[1:]:
        other_mt = hl.import_vcf(vcf, force_bgz=vcf.split('.')[-1]=='gz')
        other_mt = other_mt.semi_join_rows(common_rows)
        other_mt = other_mt.select_rows()  # drop row fields
        mt = mt.union_cols(other_mt)

    # Keep only row fields from the first VCF
    mt = mt.select_rows(*list(mt.row))

    return mt


def main():
    parser = argparse.ArgumentParser(description="Merge VCFs in Hail")
    parser.add_argument("--vcfs", required=True, help="Comma-separated list of VCF files")
    parser.add_argument("--output", required=True, help="Output path (VCF or MT). Example: merged.vcf.bgz or merged.mt")
    parser.add_argument('--mem', type=float, help='Memory (GB).')
    parser.add_argument('--genome_build', help='Genome build (GRCh37 or GRCh38).')
    
    args = parser.parse_args()

    # Split input into list
    vcf_files = args.vcfs.split(",")

    genome_build = args.genome_build
    mem = args.mem

    # Initialize Hail
    hl.init(min_block_size=128, spark_conf={
                    "spark.executor.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g"
                    }, tmp_dir="tmp", local_tmpdir="tmp",
                    default_reference=genome_build)

    # Merge
    merged_mt = merge_vcfs(vcf_files)

    first_vcf_header =hl.get_vcf_metadata(vcf_files[0])

    hl.export_vcf(merged_mt, args.output, tabix=True, metadata=first_vcf_header)

if __name__ == "__main__":
    main()
