#!/usr/bin/env python3
import hail as hl
import pandas as pd
from typing import List
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Merge vcf_files sample-wise with Hail and update INFO fields")

parser.add_argument("--vcf-files", nargs="+", help="List of VCF files to merge")
parser.add_argument("--batch-size", type=int, help="Batch size (number of vcf_files) for processing")
parser.add_argument("--output-vcf-file", type=str, help="Path for final VCF output")
parser.add_argument("--mem", type=float, help="Memory in GB to allocate for Hail processing")
parser.add_argument("--max-info-fields", nargs="+", help="List of INFO field to merge by taking the max")
parser.add_argument("--min-info-fields", nargs="+", help="List of INFO field to merge by taking the min")
parser.add_argument("--sum-info-fields", nargs="+", help="List of INFO field to merge by taking the sum")
args = parser.parse_args()

vcf_files = args.vcf_files
batch_size = args.batch_size
output_vcf = args.output_vcf_file
mem = args.mem
max_info_fields = args.max_info_fields
min_info_fields = args.min_info_fields
sum_info_fields = args.sum_info_fields

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

def import_vcf_with_partitions(
    vcf_path: str,
    partitions=None,
    reference_genome: str = "GRCh38",
) -> hl.MatrixTable:
    print(f"  Importing VCF: {vcf_path}", flush=True)
    mt = hl.import_vcf(
        vcf_path,
        reference_genome=reference_genome,
        array_elements_required=False,
        force_bgz=True,
        min_partitions=None if partitions is None else partitions,
        n_partitions=partitions,
    )
    if 'VQSLOD' in list(mt.info):
        # VQSLOD imported as str for some reason --> convert to float
        mt = mt.annotate_rows(info=mt.info.annotate(VQSLOD=hl.float(mt.info.VQSLOD)))
    return mt

def import_batch_and_union_cols(
    vcf_paths: List[str],
    reference_genome: str = "GRCh38",
) -> hl.MatrixTable:
    """
    Import a batch of vcf_files and union them by columns.
    Partitioning is derived from the first VCF only.
    """

    # Import first VCF normally
    mt = import_vcf_with_partitions(vcf_paths[0], None, reference_genome)

    # Derive partitioning from first MT
    partitions = mt.n_partitions()
    print(f"  Using {partitions} partitions for remaining imports", flush=True)

   # Import remaining vcf_files using that partitioning
    for idx, path in enumerate(vcf_paths[1:], start=2):
        print(f"  Unioning VCF {idx}/{len(vcf_paths)}")
        other = import_vcf_with_partitions(path, partitions, reference_genome)
        mt = mt.union_cols(
            other,
            row_join_type="outer",
            drop_right_row_fields=False,
        )

    print("Finished batch column union")
    return mt

def union_vcf_files_in_batches(
    vcf_paths: List[str],
    batch_size: int,
    tmp_prefix: str,
) -> List[str]:

    intermediate = []

    total_batches = (len(vcf_paths) + batch_size - 1) // batch_size
    print(f"Unioning {len(vcf_paths)} vcf_files in {total_batches} batches")

    for i in range(0, len(vcf_paths), batch_size):
        batch_idx = i // batch_size
        batch = vcf_paths[i:i + batch_size]

        print(f"\nProcessing batch {batch_idx + 1}/{total_batches}")
        merged = import_batch_and_union_cols(batch)

        out = f"{tmp_prefix}.batch_{batch_idx}.mt"
        print(f"Writing batch MT to {out}")
        merged.write(out, overwrite=True)

        intermediate.append(out)

    print("\nFinished all batch unions")
    return intermediate

def recursive_union_cols(
    mt_paths: List[str],
    final_out: str,
):
    current = mt_paths
    round_num = 0

    print(f"\nStarting recursive column union with {len(current)} MTs")

    while len(current) > 1:
        print(f"\nUnion round {round_num}: {len(current)} MTs")
        next_round = []

        for i in range(0, len(current), 2):
            if i + 1 == len(current):
                print(f"  Carrying forward MT {current[i]}")
                next_round.append(current[i])
                continue

            print(f"  Unioning MT pair {i // 2 + 1}")
            mt1 = hl.read_matrix_table(current[i])
            mt2 = hl.read_matrix_table(current[i + 1])

            merged = mt1.union_cols(
                mt2,
                row_join_type="outer",
                drop_right_row_fields=False,
            )

            out = f"{final_out}.round_{round_num}_{i // 2}.mt"
            print(f"  Writing intermediate MT to {out}")
            merged.write(out, overwrite=True)

            next_round.append(out)

        current = next_round
        round_num += 1

    print(f"\nWriting final unioned MT to {final_out}")
    hl.read_matrix_table(current[0]).write(final_out, overwrite=True)

def union_many_vcf_files(
    vcf_paths: List[str],
    batch_size: int,
    tmp_prefix: str,
    final_out: str,
):
    print("=== Starting VCF column union workflow ===")

    intermediates = union_vcf_files_in_batches(
        vcf_paths=vcf_paths,
        batch_size=batch_size,
        tmp_prefix=tmp_prefix,
    )

    recursive_union_cols(
        mt_paths=intermediates,
        final_out=final_out,
    )

    print("=== VCF column union workflow complete ===")

tmp_prefix = output_vcf.split('.vcf')[0]
final_out = f"{tmp_prefix}.merged.mt"
union_many_vcf_files(
    vcf_paths=vcf_files,
    batch_size=batch_size,
    tmp_prefix=tmp_prefix,
    final_out=final_out,
)

final_mt = hl.read_matrix_table(final_out)

# Don't care as much about the other fields --> can just leave them as the first VCF's values (default behavior)
# Need these fields updated for hard filtering and PURF model
merge_strategy = {"max": max_info_fields, "min": min_info_fields, "sum": sum_info_fields}

# Get all row fields
row_fields = list(final_mt.row)

# Filter for per-VCF INFO fields using "info_{n}" pattern (except first will just be "info")
info_fields_list = [f for f in row_fields if f.startswith("info")]
agg_map = {"min": hl.min, "max": hl.max, "sum": hl.sum}

for strat, fields in merge_strategy.items():
    agg_fn = agg_map[strat]
    
    for field in fields:       
        # aggregate into single INFO field
        final_mt = final_mt.annotate_rows(
            info = final_mt.info.annotate(**{
                field: agg_fn(hl.array([final_mt[f][field] for f in info_fields_list]))
            })
        )

## Need to update INFO field values after merging vcf_files across samples ##
# Recalculate AC, AF, AN after merge
final_mt = hl.variant_qc(final_mt)
final_mt = final_mt.annotate_rows(info=final_mt.info.annotate(
    cohort_AN=final_mt.variant_qc.AN,
    cohort_AC=final_mt.variant_qc.AC[1],
    cohort_AF=final_mt.variant_qc.AF[1]))

# Assumes uniform header across vcf_files --> just grab first VCF's header
header = hl.get_vcf_metadata(vcf_files[0])

hl.export_vcf(final_mt, output_vcf,
              metadata=header)