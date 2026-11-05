#!/usr/bin/env python3
import hail as hl
import pandas as pd
from typing import Dict, List, Optional, Tuple
import argparse
import numpy as np
from hail_vcf_utils import split_multi_ssc, subset_allele_info_fields

FILL_MISSING_MODES = {
    "all": "every missing GT becomes 0/0, including no-calls (./.) present in the input VCFs",
    "absent": "0/0 only where the sample's VCF has no record for the variant",
    "absent_at_genotyped_loci": "as 'absent', but only where the sample's VCF has some record at the same locus",
    "none": "missing GTs stay missing",
}

def parse_args():
    parser = argparse.ArgumentParser(description="Merge vcf_files sample-wise with Hail and update INFO fields")

    parser.add_argument("--vcf-files", nargs="+", help="List of VCF files to merge")
    parser.add_argument("--batch-size", type=int, default=10,
                        help="Number of VCFs unioned together before each intermediate MT is written; "
                             "the batch MTs are then unioned pairwise")
    parser.add_argument("--output-vcf-file", type=str, help="Path for final VCF output")
    parser.add_argument("--mem", type=float, help="Memory in GB to allocate for Hail processing")
    parser.add_argument("--max-info-fields", nargs="*", default=[], help="List of INFO fields to merge by taking the max")
    parser.add_argument("--min-info-fields", nargs="*", default=[], help="List of INFO fields to merge by taking the min")
    parser.add_argument("--sum-info-fields", nargs="*", default=[], help="List of INFO fields to merge by taking the sum")
    parser.add_argument("--genome-build", type=str, default="GRCh38", help="Reference genome")
    parser.add_argument("--sample-file", type=str, default=None,
                        help="OPTIONAL: file with one sample ID per line (first tab-separated column); only these samples are kept. "
                             "IDs are matched before --cohort-prefixes is applied")
    parser.add_argument("--cohort-prefixes", nargs="*", default=[],
                        help="OPTIONAL: one prefix per VCF, in --vcf-files order; samples are renamed to <prefix>:<sample>")
    parser.add_argument("--format-fields-to-drop", nargs="*", default=[],
                        help="OPTIONAL: FORMAT fields to drop (GT cannot be dropped); default keeps every field shared by all VCFs")
    parser.add_argument("--site-join", choices=["outer", "inner"], default="outer",
                        help="outer: keep variants in any VCF; inner: keep only variants in every VCF")
    parser.add_argument("--fill-missing", choices=list(FILL_MISSING_MODES), default="absent",
                        help="; ".join(f"{mode}: {desc}" for mode, desc in FILL_MISSING_MODES.items()))
    parser.add_argument("--header-file", type=str, default=None,
                        help="OPTIONAL: VCF header to use for every input VCF (e.g. when tags are missing from the input headers)")
    return parser.parse_args()

def import_vcf_raw(
    vcf_path: str,
    reference_genome: str,
    header_file: Optional[str] = None,
    partitions=None,
) -> hl.MatrixTable:
    mt = hl.import_vcf(
        path=vcf_path,
        reference_genome=reference_genome,
        header_file=header_file,
        array_elements_required=False,
        force_bgz=True,
        min_partitions=None if partitions is None else partitions
    )
    # Convert VQSLOD to float if imported as string
    if 'VQSLOD' in list(mt.info):
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                VQSLOD=hl.float(mt.info.VQSLOD)
            )
        )
    return mt

def prepare_vcf(
    vcf_input: dict,
    entry_fields: List[str],
    samples_to_keep: Optional[List[str]],
) -> hl.MatrixTable:
    """
    Prepare a VCF imported by import_and_check_vcfs() for the union: subset samples, split
    multiallelics, subset per-allele INFO, select FORMAT fields and rename samples.
    """
    print(f"  Preparing VCF: {vcf_input['path']}", flush=True)
    mt = vcf_input['mt']
    if samples_to_keep is not None:
        mt = mt.filter_cols(hl.literal(set(samples_to_keep)).contains(mt.s))
    # Split multiallelics before merging, so the same allele is keyed the same way in every VCF
    mt = split_multi_ssc(mt)
    mt = subset_allele_info_fields(mt, vcf_input['info_header'])
    # _present marks entries from this VCF, so a no-call (./.) can be told apart from a variant the VCF never had
    mt = mt.select_entries(*entry_fields, _present=True)
    # Add prefix to sample names if optioned
    if vcf_input['prefix']:
        mt = mt.key_cols_by(s=hl.str(f"{vcf_input['prefix']}:") + mt.s)
    # _vcf_idx records which input VCF each sample came from (used by fill_missing_genotypes)
    mt = mt.select_cols(_vcf_idx=vcf_input['vcf_idx'])
    # Drop fields introduced by multiallelic split
    return mt.drop('num_alleles', 'a_index', 'was_split')

def write_genotyped_loci(mt: hl.MatrixTable, out: str):
    loci = mt.rows().key_by('locus').select().distinct()
    loci.write(out, overwrite=True)

def prepare_batch_and_union_cols(
    vcf_inputs: List[dict],
    prepare_args: dict,
    row_join_type: str,
    loci_prefix: Optional[str],
) -> hl.MatrixTable:
    """
    Prepare a batch of imported VCFs and union them by columns, keeping all non-key row fields.
    Also, write out genotyped loci in each VCF if optioned by loci_prefix.
    """

    mt = prepare_vcf(vcf_inputs[0], **prepare_args)
    if loci_prefix:
        write_genotyped_loci(mt, f"{loci_prefix}_{vcf_inputs[0]['vcf_idx']}.ht")

    for idx, vcf_input in enumerate(vcf_inputs[1:], start=2):
        print(f"  Unioning VCF {idx}/{len(vcf_inputs)}")
        other = prepare_vcf(vcf_input, **prepare_args)
        if loci_prefix:
            write_genotyped_loci(other, f"{loci_prefix}_{vcf_input['vcf_idx']}.ht")
        mt = mt.union_cols(
            other,
            row_join_type=row_join_type,
            drop_right_row_fields=False,
        )

    print("Finished batch column union")
    return mt

def union_vcf_files_in_batches(
    vcf_inputs: List[dict],
    batch_size: int,
    tmp_prefix: str,
    prepare_args: dict,
    row_join_type: str,
    loci_prefix: Optional[str],
) -> List[str]:

    intermediate = []

    total_batches = (len(vcf_inputs) + batch_size - 1) // batch_size
    print(f"Unioning {len(vcf_inputs)} VCFs in {total_batches} batches")

    for i in range(0, len(vcf_inputs), batch_size):
        batch_idx = i // batch_size
        batch = vcf_inputs[i:i + batch_size]

        print(f"\nProcessing batch {batch_idx + 1}/{total_batches}")
        merged = prepare_batch_and_union_cols(batch, prepare_args, row_join_type, loci_prefix)

        out = f"{tmp_prefix}.batch_{batch_idx}.mt"
        print(f"Writing batch MT to {out}")
        merged.write(
            out,
            overwrite=True
        )

        intermediate.append(out)

    print("\nFinished all batch unions")
    return intermediate

def recursive_union_cols(
    mt_paths: List[str],
    tmp_prefix: str,
    row_join_type: str,
) -> str:
    """
    Union the batch MTs pairwise, writing each round's unions, until one MT is left.
    Returns the path of that MT (a batch MT if there is only one batch).
    """
    current = mt_paths
    round_num = 0

    print(f"\nStarting recursive column union with {len(current)} MTs")

    while len(current) > 1:
        print(f"\nUnion round {round_num}: {len(current)} MTs")
        next_round = []

        for i in range(0, len(current), 2):
            if i + 1 == len(current):
                # Odd MT left over, carry forward
                print(f"  Carrying forward MT {current[i]}")
                next_round.append(current[i])
                continue

            # Union MT pairs
            print(f"  Unioning MT pair {i // 2 + 1}")
            mt1 = hl.read_matrix_table(path=current[i])
            mt2 = hl.read_matrix_table(path=current[i + 1])

            merged = mt1.union_cols(
                mt2,
                row_join_type=row_join_type,
                drop_right_row_fields=False,
            )

            out = f"{tmp_prefix}.round_{round_num}_{i // 2}.mt"
            print(f"  Writing intermediate MT to {out}")
            merged.write(
                out,
                overwrite=True
            )

            next_round.append(out)

        current = next_round
        round_num += 1

    # The last MT written already holds the full union, so it is used as is rather than copied
    print(f"\nFinal unioned MT: {current[0]}")
    return current[0]

def union_many_vcf_files(
    vcf_inputs: List[dict],
    batch_size: int,
    tmp_prefix: str,
    prepare_args: dict,
    row_join_type: str,
    loci_prefix: Optional[str],
) -> str:
    """
    Union all input VCFs by columns; returns the path of the unioned MT.
    """
    print("=== Starting VCF column union workflow ===")

    intermediates = union_vcf_files_in_batches(
        vcf_inputs=vcf_inputs,
        batch_size=batch_size,
        tmp_prefix=tmp_prefix,
        prepare_args=prepare_args,
        row_join_type=row_join_type,
        loci_prefix=loci_prefix,
    )

    final_path = recursive_union_cols(
        mt_paths=intermediates,
        tmp_prefix=tmp_prefix,
        row_join_type=row_join_type,
    )

    print("=== VCF column union workflow complete ===")
    return final_path

def fill_missing_genotypes(
    mt: hl.MatrixTable,
    fill_missing: str,
    loci_ht_paths: Dict[int, str],
) -> hl.MatrixTable:
    hom_ref = hl.call(0, 0)
    if fill_missing == "all":
        mt = mt.annotate_entries(GT=hl.coalesce(mt.GT, hom_ref))
    elif fill_missing == "absent":
        mt = mt.annotate_entries(GT=hl.if_else(hl.is_defined(mt._present), mt.GT, hom_ref))
    elif fill_missing == "absent_at_genotyped_loci":
        # Arrays of inputs (by _vcf_idx) with any record at each locus in the MT
        genotyped_vcf_idxs = hl.array([
            hl.or_missing(hl.is_defined(hl.read_table(path)[mt.locus]), vcf_idx)
            for vcf_idx, path in sorted(loci_ht_paths.items())
        ])
        mt = mt.annotate_rows(_genotyped_vcf_idxs=hl.set(genotyped_vcf_idxs.filter(hl.is_defined)))
        # If sample-variant entry is in the input VCF for the sample, use the GT from
        # the input VCF; else if the locus was genotyped but there was no record for
        # this variant, fill with 0/0
        mt = mt.annotate_entries(GT=hl.if_else(
            hl.is_defined(mt._present),
            mt.GT,
            hl.or_missing(mt._genotyped_vcf_idxs.contains(mt._vcf_idx), hom_ref)
        ))
        mt = mt.drop('_genotyped_vcf_idxs')
    return mt

def read_sample_file(sample_file: str) -> List[str]:
    with hl.hadoop_open(sample_file, 'r') as f:
        return [line.strip().split('\t')[0] for line in f if line.strip()]

def import_and_check_vcfs(
    vcf_paths: List[str],
    reference_genome: str,
    header_file: Optional[str],
    samples_to_keep: Optional[List[str]],
    cohort_prefixes: List[str],
    format_fields_to_drop: List[str],
) -> Tuple[List[dict], List[str]]:
    """
    Import every VCF once (lazily; partitioning is derived from the first VCF) and check headers
    and samples before any data is processed. Returns the inputs to merge, each with its imported
    MatrixTable and INFO header (VCFs left with no samples after subsetting are skipped), and the
    FORMAT fields to keep: those present in every VCF with the same type and not in
    format_fields_to_drop, in the order of the first VCF with samples to keep.
    Fails on duplicate sample IDs, so they never reach union_cols.
    """
    if cohort_prefixes and len(cohort_prefixes) != len(vcf_paths):
        raise ValueError(f"--cohort-prefixes has {len(cohort_prefixes)} values for {len(vcf_paths)} VCFs")
    if 'GT' in format_fields_to_drop:
        raise ValueError("GT cannot be dropped: it is needed to split multiallelics and recompute AC/AF/AN")

    # Samples in input sample file to subset VCFs if optioned, otherwise disregarded
    keep = set(samples_to_keep) if samples_to_keep is not None else None
    # Initialize VCF inputs to pass along for merging
    inputs = []
    # Initialize types of entry fields which match across all VCFs, in the order of the first VCF
    # with samples to keep (dicts keep insertion order, and narrowing below never reorders)
    entry_types = None
    # Entry fields missing from, or with a different type in, at least one VCF
    unshared_fields = set()
    # If samples_to_keep is optioned, samples_found holds those samples that were found
    # and subsetted to in the VCFs
    samples_found = set()
    # The header file's INFO metadata applies to every VCF; otherwise each VCF's own is read below
    shared_info_header = hl.get_vcf_metadata(header_file)['info'] if header_file else None
    partitions = None
    for vcf_idx, path in enumerate(vcf_paths):
        print(f"  Importing VCF: {path}", flush=True)
        mt = import_vcf_raw(path, reference_genome, header_file, partitions)
        if partitions is None:
            # Derive partitioning from first MT
            partitions = mt.n_partitions()
            print(f"  Using {partitions} partitions for remaining imports", flush=True)
        samples = mt.s.collect()
        if keep is not None:
            samples = [s for s in samples if s in keep]
            samples_found.update(samples)
            if not samples:
                print(f"  Skipping {path}: none of its samples are in the sample file", flush=True)
                continue
        prefix = cohort_prefixes[vcf_idx] if cohort_prefixes else None
        inputs.append({'path': path, 'vcf_idx': vcf_idx, 'prefix': prefix, 'mt': mt,
                       'info_header': shared_info_header or hl.get_vcf_metadata(path)['info'],
                       'samples': [f"{prefix}:{s}" if prefix else s for s in samples]})

        types = {f: mt.entry[f].dtype for f in mt.entry}
        if entry_types is None:
            entry_types = types
        else:
            unshared_fields.update(f for f in set(entry_types) | set(types) if entry_types.get(f) != types.get(f))
            entry_types = {f: t for f, t in entry_types.items() if types.get(f) == t}

    if not inputs:
        raise ValueError("No samples left to merge")
    if keep is not None:
        not_found = sorted(keep - samples_found)
        if not_found:
            print(f"  {len(not_found)} samples in the sample file are in no VCF, e.g. {not_found[:10]}", flush=True)

    all_samples = [s for inp in inputs for s in inp['samples']]
    duplicates = sorted(s for s, n in pd.Series(all_samples).value_counts().items() if n > 1)
    if duplicates:
        raise ValueError(f"{len(duplicates)} sample IDs appear in more than one VCF, e.g. {duplicates[:10]}; "
                         f"use --cohort-prefixes to make them unique")

    unknown = [f for f in format_fields_to_drop if f not in entry_types]
    if unknown:
        print(f"  FORMAT fields to drop that are not in every VCF with the same type (already dropped): {unknown}", flush=True)
    entry_fields = [f for f in entry_types if f not in set(format_fields_to_drop)]
    if 'GT' not in entry_fields:
        raise ValueError("GT is not a FORMAT field of every VCF")
    if unshared_fields:
        print(f"  Dropping FORMAT fields not shared by every VCF with the same type: {sorted(unshared_fields)}", flush=True)
    requested = [f for f in format_fields_to_drop if f in entry_types]
    if requested:
        print(f"  Dropping FORMAT fields as requested: {requested}", flush=True)

    print(f"  Merging {len(all_samples)} samples from {len(inputs)} VCFs", flush=True)
    return inputs, entry_fields

def main():
    args = parse_args()

    vcf_files = args.vcf_files
    batch_size = args.batch_size
    output_vcf = args.output_vcf_file
    mem = args.mem
    max_info_fields = args.max_info_fields
    min_info_fields = args.min_info_fields
    sum_info_fields = args.sum_info_fields
    genome_build = args.genome_build
    header_file = args.header_file
    site_join = args.site_join
    fill_missing = args.fill_missing

    hl.init(
        min_block_size=128,
        local=f"local[*]",
        spark_conf={
            "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
            "spark.speculation": 'true'
        },
        tmp_dir="tmp",
        local_tmpdir="tmp",
    )

    samples_to_keep = read_sample_file(args.sample_file) if args.sample_file else None

    # -------------------------------
    # Check headers, samples and FORMAT fields before merging
    # -------------------------------
    vcf_inputs, entry_fields = import_and_check_vcfs(
        vcf_paths=vcf_files,
        reference_genome=genome_build,
        header_file=header_file,
        samples_to_keep=samples_to_keep,
        cohort_prefixes=args.cohort_prefixes,
        format_fields_to_drop=args.format_fields_to_drop,
    )
    prepare_args = {
        'entry_fields': entry_fields,
        'samples_to_keep': samples_to_keep,
    }

    tmp_prefix = output_vcf.split('.vcf')[0]
    loci_prefix = f"{tmp_prefix}.genotyped_loci" if fill_missing == "absent_at_genotyped_loci" else None
    final_path = union_many_vcf_files(
        vcf_inputs=vcf_inputs,
        batch_size=batch_size,
        tmp_prefix=tmp_prefix,
        prepare_args=prepare_args,
        row_join_type=site_join,
        loci_prefix=loci_prefix,
    )

    final_mt = hl.read_matrix_table(path=final_path)

    # -------------------------------
    # Fill missing GTs (see FILL_MISSING_MODES)
    # -------------------------------
    print(f"Filling missing GTs: {fill_missing} ({FILL_MISSING_MODES[fill_missing]})", flush=True)
    loci_ht_paths = {inp['vcf_idx']: f"{loci_prefix}_{inp['vcf_idx']}.ht" for inp in vcf_inputs} if loci_prefix else {}
    final_mt = fill_missing_genotypes(final_mt, fill_missing, loci_ht_paths)
    # _present and _vcf_idx were only needed to fill missing GTs
    final_mt = final_mt.drop('_present', '_vcf_idx')

    # -------------------------------
    # Update INFO field values after merging VCFs across samples
    # -------------------------------
    merge_strategy = {"max": max_info_fields, "min": min_info_fields, "sum": sum_info_fields}

    # -------------------------------
    # Collect all row fields
    # -------------------------------
    row_fields = list(final_mt.row)

    # -------------------------------
    # QUAL aggregation (min across VCFs)
    # -------------------------------
    qual_fields_list = [f for f in row_fields if f.startswith("qual")]
    final_mt = final_mt.annotate_rows(
        qual=hl.min(
            hl.array([final_mt[f] for f in qual_fields_list])
        )
    )

    # -------------------------------
    # FILTER aggregation (union across VCFs)
    # -------------------------------
    filters_fields_list = [f for f in row_fields if f.startswith("filters")]
    final_mt = final_mt.annotate_rows(
        filters=hl.flatten(
            hl.set([final_mt[f] for f in filters_fields_list])
        )
    )

    # -------------------------------
    # ID aggregation (first non-missing across VCFs)
    # -------------------------------
    rsid_fields_list = [f for f in row_fields if f.startswith("rsid")]
    if rsid_fields_list:
        final_mt = final_mt.annotate_rows(
            rsid=hl.coalesce(*[final_mt[f] for f in rsid_fields_list])
        )

    # -------------------------------
    # INFO fields aggregation
    # -------------------------------
    info_fields_list = [f for f in row_fields if f.startswith("info")]

    agg_map = {"min": hl.min, "max": hl.max, "sum": hl.sum}

    # Start with empty struct to avoid issues when original 'info' is empty
    final_mt = final_mt.annotate_rows(
        info_new=hl.struct()
    )

    # Step 1: Aggregate merge_strategy fields (min/max/sum)
    for strat, fields in merge_strategy.items():
        agg_fn = agg_map[strat]
        for field in fields:
            # Only include INFO structs that actually have this field
            exprs = [
                final_mt[f][field]
                for f in info_fields_list
                if field in final_mt[f].dtype
            ]
            if not exprs:
                print(f"  Skipping INFO field {field} ({strat}): not in any VCF", flush=True)
                continue
            final_mt = final_mt.annotate_rows(
                info_new=final_mt.info_new.annotate(
                    **{
                        field: agg_fn(
                            hl.array(exprs)
                        )
                    }
                )
            )

    # Step 2: For all other INFO fields, take first non-missing value
    # (kept in input header order, not set order, so the output is the same on every run)

    all_info_keys = []
    for f in info_fields_list:
        all_info_keys += [k for k in final_mt[f].dtype.keys() if k not in all_info_keys]

    handled_keys = set(sum(merge_strategy.values(), []))
    other_info_keys = [k for k in all_info_keys if k not in handled_keys]

    for field in other_info_keys:
        # Only include INFO structs that actually have this field
        exprs = [
            final_mt[f][field]
            for f in info_fields_list
            if field in final_mt[f].dtype
        ]

        final_mt = final_mt.annotate_rows(
            info_new=final_mt.info_new.annotate(
                **{
                    field: hl.find(lambda x: hl.is_defined(x), exprs)
                }
            )
        )

    # -------------------------------
    # Checkpoint before cohort-level QC
    # -------------------------------
    final_out_updated = f"{tmp_prefix}.merged.updated.info.mt"
    final_mt = final_mt.checkpoint(final_out_updated)

    # -------------------------------
    # Recalculate cohort-level QC after filling missing GTs
    # -------------------------------
    final_mt = hl.variant_qc(final_mt)
    final_mt = final_mt.annotate_rows(
        info_new=final_mt.info_new.annotate(
            cohort_AN=final_mt.variant_qc.AN,
            cohort_AC=final_mt.variant_qc.AC[1],
            cohort_AF=final_mt.variant_qc.AF[1]
        )
    )
    # AC/AF/AN were taken from the first VCF with the site; replace them with the merged counts
    merged_counts = {
        "AC": [final_mt.variant_qc.AC[1]],
        "AF": [final_mt.variant_qc.AF[1]],
        "AN": final_mt.variant_qc.AN,
    }
    final_mt = final_mt.annotate_rows(
        info_new=final_mt.info_new.annotate(
            **{field: value for field, value in merged_counts.items() if field in final_mt.info_new}
        )
    )

    # Replace info with info_new for exporting to VCF
    final_mt = final_mt.annotate_rows(
        info=final_mt.info_new
    )
    final_mt = final_mt.select_rows(*[f for f in ('rsid', 'qual', 'filters') if f in final_mt.row], 'info')

    # Assumes uniform header across VCFs --> grab the header file, or else the first VCF's header
    header = hl.get_vcf_metadata(header_file or vcf_files[0])
    header['info'].update({
        'cohort_AC': {'Description': 'Alternate allele count over all merged samples', 'Number': '1', 'Type': 'Integer'},
        'cohort_AF': {'Description': 'Alternate allele frequency over all merged samples', 'Number': '1', 'Type': 'Float'},
        'cohort_AN': {'Description': 'Total number of called alleles over all merged samples', 'Number': '1', 'Type': 'Integer'},
    })

    hl.export_vcf(
        final_mt,
        output=output_vcf,
        metadata=header,
        tabix=output_vcf.endswith(('.bgz', '.gz'))
    )

if __name__ == "__main__":
    main()
