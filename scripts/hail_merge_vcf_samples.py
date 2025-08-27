#!/usr/bin/env python3

import argparse
import numpy as np
import hail as hl

def intersect_entry_fields(mts):
    """
    Get entry field names shared across all MatrixTables with identical types.
    """
    common_fields = set(mts[0].entry.keys())
    for mt in mts[1:]:
        common_fields &= set(mt.entry.keys())

    # Check types and keep only fields with matching types
    consistent_fields = []
    for f in common_fields:
        t0 = mts[0].entry[f].dtype
        if all(mt.entry[f].dtype == t0 for mt in mts[1:]):
            consistent_fields.append(f)

    return consistent_fields

def merge_vcfs(vcf_files):
    # Load first VCF
    mt = hl.import_vcf(vcf_files[0], force_bgz=vcf_files[0].split('.')[-1]=='gz')

    # Collect common rows across all VCFs
    mts = [mt]
    common_rows = mt.rows()
    for vcf in vcf_files[1:]:
        other_mt = hl.import_vcf(vcf, force_bgz=vcf.split('.')[-1]=='gz')
        mts.append(other_mt)
        common_rows = common_rows.semi_join(other_mt.rows())

    # Restrict all datasets to common rows
    mts = [m.semi_join_rows(common_rows) for m in mts]

    # Find entry fields in common with matching dtypes
    common_fields = set(mts[0].entry.keys())
    for m in mts[1:]:
        common_fields &= set(m.entry.keys())

    consistent_fields = []
    for f in common_fields:
        dtype0 = mts[0].entry[f].dtype
        if all(m.entry[f].dtype == dtype0 for m in mts[1:]):
            consistent_fields.append(f)

    # Keep only consistent entry fields
    mts = [m.select_entries(*consistent_fields) for m in mts]

    # Merge while keeping only row fields from the first VCF
    merged = mts[0]
    for m in mts[1:]:
        m = m.select_rows()  # drop row fields
        merged = merged.union_cols(m)

    # Keep only the non-key row fields from the first VCF
    row_keys = list(merged.row_key)  # ['locus', 'alleles']
    row_fields = [f for f in mts[0].row if f not in row_keys]
    merged = merged.select_rows(*row_fields)

    return merged

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
