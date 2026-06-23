version 1.0

workflow VcfToBed {

  input {
    File    vcf_file
    String  cohort_prefix
    String  variant_interpretation_docker
    File    outlier_samples_file
    File    ped_file

    Array[String] svtypes            = ["DEL", "DUP", "CNV"]
    Float   gnomad_af_threshold      = 0.01
    Int     min_proband_ac           = 1
    Int     max_proband_ac           = 5
    Boolean filter_predicted_lof                 = true
    Boolean filter_predicted_intragenic_exon_dup = true
    Boolean filter_predicted_copy_gain           = true
    Boolean filter_predicted_partial_exon_dup    = true
  }

  call FilterAndSplitVcf {
    input:
      vcf_file                             = vcf_file,
      cohort_prefix                        = cohort_prefix,
      variant_interpretation_docker        = variant_interpretation_docker,
      outlier_samples_file                 = outlier_samples_file,
      ped_file                             = ped_file,
      svtypes                              = svtypes,
      gnomad_af_threshold                  = gnomad_af_threshold,
      min_proband_ac                       = min_proband_ac,
      max_proband_ac                       = max_proband_ac,
      filter_predicted_lof                 = filter_predicted_lof,
      filter_predicted_intragenic_exon_dup = filter_predicted_intragenic_exon_dup,
      filter_predicted_copy_gain           = filter_predicted_copy_gain,
      filter_predicted_partial_exon_dup    = filter_predicted_partial_exon_dup
  }

  output {
    # Outlier-exclusive outputs (main PROBAND_AC branch)
    File outlier_vcf                              = FilterAndSplitVcf.outlier_vcf
    File outlier_vcf_tbi                          = FilterAndSplitVcf.outlier_vcf_tbi
    File outlier_bed                              = FilterAndSplitVcf.outlier_bed

    # Non-outlier outputs (main PROBAND_AC branch, outlier samples fully excluded)
    File nonoutlier_vcf                           = FilterAndSplitVcf.nonoutlier_vcf
    File nonoutlier_vcf_tbi                       = FilterAndSplitVcf.nonoutlier_vcf_tbi
    File nonoutlier_bed                           = FilterAndSplitVcf.nonoutlier_bed

    # All-samples BED (main PROBAND_AC branch)
    File all_samples_bed                          = FilterAndSplitVcf.all_samples_bed

    # Variant ID lists (main branch)
    File ids_outlier_exclusive                    = FilterAndSplitVcf.ids_outlier_exclusive
    File ids_in_nonoutliers                       = FilterAndSplitVcf.ids_in_nonoutliers

    # Parent-only outputs: PROBAND_AC=0 branch, split by outlier status
    File no_proband_carriers_outlier_vcf          = FilterAndSplitVcf.no_proband_carriers_outlier_vcf
    File no_proband_carriers_outlier_vcf_tbi      = FilterAndSplitVcf.no_proband_carriers_outlier_vcf_tbi
    File no_proband_carriers_outlier_bed          = FilterAndSplitVcf.no_proband_carriers_outlier_bed
    File no_proband_carriers_nonoutlier_vcf       = FilterAndSplitVcf.no_proband_carriers_nonoutlier_vcf
    File no_proband_carriers_nonoutlier_vcf_tbi   = FilterAndSplitVcf.no_proband_carriers_nonoutlier_vcf_tbi
    File no_proband_carriers_nonoutlier_bed       = FilterAndSplitVcf.no_proband_carriers_nonoutlier_bed
    File ids_no_proband_carriers_outlier          = FilterAndSplitVcf.ids_no_proband_carriers_outlier
    File ids_no_proband_carriers_nonoutlier       = FilterAndSplitVcf.ids_no_proband_carriers_nonoutlier
    File proband_ids                              = FilterAndSplitVcf.proband_ids
  }
}

task FilterAndSplitVcf {

  input {
    File          vcf_file
    String        cohort_prefix
    String        variant_interpretation_docker
    File          outlier_samples_file
    File          ped_file
    Array[String] svtypes                              = ["DEL", "DUP", "CNV"]
    Float         gnomad_af_threshold                  = 0.01
    Int           min_proband_ac                       = 1
    Int           max_proband_ac                       = 5
    Boolean       filter_predicted_lof                 = true
    Boolean       filter_predicted_intragenic_exon_dup = true
    Boolean       filter_predicted_copy_gain           = true
    Boolean       filter_predicted_partial_exon_dup    = true
    Int           mem_gb                               = 8
    Int           cpu_cores                            = 2
    Int           preemptible_tries                    = 3
    Int           max_retries                          = 1
    Int           boot_disk_gb                         = 10
  }

  Float input_size = size(vcf_file, "GB")
  Int   disk_gb    = ceil(10.0 + input_size * 10.0)

  runtime {
    memory:         "~{mem_gb} GB"
    disks:          "local-disk ~{disk_gb} HDD"
    cpu:            cpu_cores
    preemptible:    preemptible_tries
    maxRetries:     max_retries
    docker:         variant_interpretation_docker
    bootDiskSizeGb: boot_disk_gb
  }

  command <<<
    set -eou pipefail

    VCF="~{vcf_file}"
    OUTLIERS="~{outlier_samples_file}"
    PED="~{ped_file}"
    PREFIX="~{cohort_prefix}"

    # ---------------------------------------------------------------
    # Step 0: Extract proband IDs from PED file
    #         PED format: family_id, sample_id, father_id, mother_id,
    #                     sex, affected (2=affected/proband)
    #         Primary:  affected status = 2
    #         Fallback: samples with both parents listed (col3/col4 != 0)
    # ---------------------------------------------------------------

    awk '$6 == 2 {print $2}' "$PED" | sort -u > proband_ids.txt

    if [ ! -s proband_ids.txt ]; then
      echo "WARNING: No affected=2 samples found in PED; falling back to samples with both parents listed"
      awk '$3 != "0" && $3 != "" && $4 != "0" && $4 != "" {print $2}' "$PED" | sort -u > proband_ids.txt
    fi

    echo "Probands identified from PED: $(wc -l < proband_ids.txt)"

    # ---------------------------------------------------------------
    # Step 1: Restrict to requested SVTYPEs (all FILTER statuses)
    # ---------------------------------------------------------------
    SVTYPE_ARRAY=(~{sep=" " svtypes})

    SVTYPE_CLAUSES=()
    for ST in "${SVTYPE_ARRAY[@]}"; do
      SVTYPE_CLAUSES+=("SVTYPE=\"${ST}\"")
    done
    SVTYPE_EXPR=$(printf '%s | ' "${SVTYPE_CLAUSES[@]}")
    SVTYPE_EXPR="(${SVTYPE_EXPR% | })"

    bcftools view \
      -f "" \
      -i "${SVTYPE_EXPR}" \
      "$VCF" \
      -Oz -o svtype_filtered.vcf.gz

    bcftools index -t svtype_filtered.vcf.gz

    # ---------------------------------------------------------------
    # Step 2: Build shared functional + AF expression
    #         PROBAND_AC thresholds are applied separately per branch
    # ---------------------------------------------------------------
    FUNC_CLAUSES=()

    if [ "~{filter_predicted_lof}" = "true" ]; then
      FUNC_CLAUSES+=('INFO/PREDICTED_LOF!="."')
    fi
    if [ "~{filter_predicted_intragenic_exon_dup}" = "true" ]; then
      FUNC_CLAUSES+=('INFO/PREDICTED_INTRAGENIC_EXON_DUP!="."')
    fi
    if [ "~{filter_predicted_copy_gain}" = "true" ]; then
      FUNC_CLAUSES+=('INFO/PREDICTED_COPY_GAIN!="."')
    fi
    if [ "~{filter_predicted_partial_exon_dup}" = "true" ]; then
      FUNC_CLAUSES+=('INFO/PREDICTED_PARTIAL_EXON_DUP!="."')
    fi

    if [ ${#FUNC_CLAUSES[@]} -eq 0 ]; then
      FUNC_EXPR="1=1"
    else
      FUNC_EXPR=$(printf '%s | ' "${FUNC_CLAUSES[@]}")
      FUNC_EXPR="(${FUNC_EXPR% | })"
    fi

    AF_EXPR="INFO/gnomad_v4.1_sv_AF<~{gnomad_af_threshold}"

    # Main branch: variants with proband carriers (PROBAND_AC min-max)
    FILTER_EXPR="${FUNC_EXPR} & INFO/PROBAND_AC>=~{min_proband_ac} & INFO/PROBAND_AC<=~{max_proband_ac} & ${AF_EXPR}"

    # Parent-only branch: variants with zero proband carriers (PROBAND_AC=0)
    PARENT_FILTER_EXPR="${FUNC_EXPR} & INFO/PROBAND_AC=0 & ${AF_EXPR}"

    echo "SVTYPE filter:         ${SVTYPE_EXPR}"
    echo "Main variant filter:   ${FILTER_EXPR}"
    echo "Parent variant filter: ${PARENT_FILTER_EXPR}"

    # ---------------------------------------------------------------
    # Step 3a: Main branch — apply PROBAND_AC min-max filter
    # ---------------------------------------------------------------
    bcftools view \
      -f "" \
      -i "${FILTER_EXPR}" \
      svtype_filtered.vcf.gz \
      -Oz -o filtered.vcf.gz

    bcftools index -t filtered.vcf.gz

    echo "Main branch variants after filtering: $(bcftools view -H filtered.vcf.gz | wc -l)"

    # ---------------------------------------------------------------
    # Step 3b: Parent-only branch — apply PROBAND_AC=0 filter
    # ---------------------------------------------------------------
    bcftools view \
      -f "" \
      -i "${PARENT_FILTER_EXPR}" \
      svtype_filtered.vcf.gz \
      -Oz -o filtered_parent_only.vcf.gz

    bcftools index -t filtered_parent_only.vcf.gz

    echo "Parent-only branch variants after filtering: $(bcftools view -H filtered_parent_only.vcf.gz | wc -l)"

    # ---------------------------------------------------------------
    # Step 4: Main branch — determine outlier-exclusive vs non-outlier
    #         variant IDs from filtered.vcf.gz
    # ---------------------------------------------------------------

    # Variant IDs carried by at least one outlier sample
    bcftools view -f "" -S "$OUTLIERS" filtered.vcf.gz \
      | bcftools view -f "" -i 'GT="alt"' \
      | bcftools query -f '%ID\n' \
      | sort -u > ids_in_outliers.txt

    # Variant IDs carried by at least one non-outlier sample
    bcftools view -f "" -S ^"$OUTLIERS" filtered.vcf.gz \
      | bcftools view -f "" -i 'GT="alt"' \
      | bcftools query -f '%ID\n' \
      | sort -u > ids_in_nonoutliers.txt

    # Outlier-exclusive: in outliers but zero non-outlier carriers
    comm -23 ids_in_outliers.txt ids_in_nonoutliers.txt > ids_outlier_exclusive.txt

    echo "Variants carried by outlier samples:     $(wc -l < ids_in_outliers.txt)"
    echo "Variants carried by non-outlier samples: $(wc -l < ids_in_nonoutliers.txt)"
    echo "Variants exclusive to outlier samples:   $(wc -l < ids_outlier_exclusive.txt)"

    # ---------------------------------------------------------------
    # Step 5: Main branch — generate outlier-exclusive and non-outlier VCFs
    # ---------------------------------------------------------------

    # Outlier-exclusive VCF — variants with zero non-outlier carriers
    bcftools view -f "" \
      -i "ID=@ids_outlier_exclusive.txt" \
      -S "$OUTLIERS" \
      filtered.vcf.gz \
      -Oz -o "${PREFIX}.outlier_samples.vcf.gz"
    bcftools index -t "${PREFIX}.outlier_samples.vcf.gz"

    # Non-outlier VCF — variants carried by at least one non-outlier,
    # samples column restricted to non-outliers only
    bcftools view -f "" \
      -i "ID=@ids_in_nonoutliers.txt" \
      -S ^"$OUTLIERS" \
      filtered.vcf.gz \
      -Oz -o "${PREFIX}.nonoutlier_samples.vcf.gz"
    bcftools index -t "${PREFIX}.nonoutlier_samples.vcf.gz"

    # ---------------------------------------------------------------
    # Step 6: Parent-only branch — split by outlier status
    #         Uses filtered_parent_only.vcf.gz (PROBAND_AC=0)
    # ---------------------------------------------------------------

    # Parent-only variant IDs carried by at least one outlier sample
    bcftools view -f "" -S "$OUTLIERS" filtered_parent_only.vcf.gz \
      | bcftools view -f "" -i 'GT="alt"' \
      | bcftools query -f '%ID\n' \
      | sort -u > ids_parent_only_in_outliers.txt

    # Parent-only variant IDs carried by at least one non-outlier sample
    bcftools view -f "" -S ^"$OUTLIERS" filtered_parent_only.vcf.gz \
      | bcftools view -f "" -i 'GT="alt"' \
      | bcftools query -f '%ID\n' \
      | sort -u > ids_parent_only_in_nonoutliers.txt

    echo "Parent-only variants in outlier samples:     $(wc -l < ids_parent_only_in_outliers.txt)"
    echo "Parent-only variants in non-outlier samples: $(wc -l < ids_parent_only_in_nonoutliers.txt)"

    # Outlier parent-only VCF — samples restricted to outliers only
    bcftools view -f "" \
      -i "ID=@ids_parent_only_in_outliers.txt" \
      -S "$OUTLIERS" \
      filtered_parent_only.vcf.gz \
      -Oz -o "${PREFIX}.no_proband_carriers_outlier.vcf.gz"
    bcftools index -t "${PREFIX}.no_proband_carriers_outlier.vcf.gz"

    # Nonoutlier parent-only VCF — samples restricted to non-outliers only
    bcftools view -f "" \
      -i "ID=@ids_parent_only_in_nonoutliers.txt" \
      -S ^"$OUTLIERS" \
      filtered_parent_only.vcf.gz \
      -Oz -o "${PREFIX}.no_proband_carriers_nonoutlier.vcf.gz"
    bcftools index -t "${PREFIX}.no_proband_carriers_nonoutlier.vcf.gz"

    # ---------------------------------------------------------------
    # Step 7: Convert all VCFs to BED
    # ---------------------------------------------------------------
    svtk vcf2bed --include-filters -i ALL \
      filtered.vcf.gz "${PREFIX}.all_samples.bed.gz"

    svtk vcf2bed --include-filters -i ALL \
      "${PREFIX}.outlier_samples.vcf.gz" "${PREFIX}.outlier_samples.bed.gz"

    svtk vcf2bed --include-filters -i ALL \
      "${PREFIX}.nonoutlier_samples.vcf.gz" "${PREFIX}.nonoutlier_samples.bed.gz"

    svtk vcf2bed --include-filters -i ALL \
      "${PREFIX}.no_proband_carriers_outlier.vcf.gz" "${PREFIX}.no_proband_carriers_outlier.bed.gz"

    svtk vcf2bed --include-filters -i ALL \
      "${PREFIX}.no_proband_carriers_nonoutlier.vcf.gz" "${PREFIX}.no_proband_carriers_nonoutlier.bed.gz"

    # ---------------------------------------------------------------
    # Step 8: Build per-subset carrier maps and fix samples column in BEDs
    #         - Main outlier BED:               outlier carriers only
    #         - Main nonoutlier BED:            non-outlier carriers only
    #         - Parent-only outlier BED:        outlier carriers only (from parent-only VCF)
    #         - Parent-only nonoutlier BED:     non-outlier carriers only (from parent-only VCF)
    # ---------------------------------------------------------------

    # Main branch: outlier-only carrier map
    bcftools view -f "" -S "$OUTLIERS" filtered.vcf.gz \
    | bcftools query -f '[%ID\t%SAMPLE\n]' -i 'GT="alt"' \
    | awk '{carriers[$1] = (carriers[$1] == "" ? $2 : carriers[$1] "," $2)}
           END {for (id in carriers) print id"\t"carriers[id]}' \
    | sort > outlier_carrier_map.txt

    echo "Main outlier carrier map: $(wc -l < outlier_carrier_map.txt) variants"

    # Main branch: non-outlier-only carrier map
    bcftools view -f "" -S ^"$OUTLIERS" filtered.vcf.gz \
    | bcftools query -f '[%ID\t%SAMPLE\n]' -i 'GT="alt"' \
    | awk '{carriers[$1] = (carriers[$1] == "" ? $2 : carriers[$1] "," $2)}
           END {for (id in carriers) print id"\t"carriers[id]}' \
    | sort > nonoutlier_carrier_map.txt

    echo "Main nonoutlier carrier map: $(wc -l < nonoutlier_carrier_map.txt) variants"

    # Parent-only branch: outlier-only carrier map
    bcftools view -f "" -S "$OUTLIERS" filtered_parent_only.vcf.gz \
    | bcftools query -f '[%ID\t%SAMPLE\n]' -i 'GT="alt"' \
    | awk '{carriers[$1] = (carriers[$1] == "" ? $2 : carriers[$1] "," $2)}
           END {for (id in carriers) print id"\t"carriers[id]}' \
    | sort > parent_outlier_carrier_map.txt

    echo "Parent-only outlier carrier map: $(wc -l < parent_outlier_carrier_map.txt) variants"

    # Parent-only branch: non-outlier-only carrier map
    bcftools view -f "" -S ^"$OUTLIERS" filtered_parent_only.vcf.gz \
    | bcftools query -f '[%ID\t%SAMPLE\n]' -i 'GT="alt"' \
    | awk '{carriers[$1] = (carriers[$1] == "" ? $2 : carriers[$1] "," $2)}
           END {for (id in carriers) print id"\t"carriers[id]}' \
    | sort > parent_nonoutlier_carrier_map.txt

    echo "Parent-only nonoutlier carrier map: $(wc -l < parent_nonoutlier_carrier_map.txt) variants"

    # Fix BED samples columns using the appropriate carrier map per BED
    python3 << 'PYEOF'
import gzip
import csv

def load_carrier_map(path):
    carrier_map = {}
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                carrier_map[parts[0]] = parts[1]
    print(f"Loaded carrier map from {path}: {len(carrier_map)} variants")
    return carrier_map

def fix_bed(bed_in, bed_out, carrier_map):
    n_fixed = 0
    n_total = 0
    with gzip.open(bed_in, 'rt') as fin, gzip.open(bed_out, 'wt') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t', lineterminator='\n')
        for i, row in enumerate(reader):
            if i == 0:
                writer.writerow(row)
                continue
            n_total += 1
            variant_id = row[3]  # name column (0-indexed)
            if variant_id in carrier_map:
                row[5] = carrier_map[variant_id]  # replace samples column
                n_fixed += 1
            writer.writerow(row)
    print(f"{bed_in}: fixed {n_fixed}/{n_total} rows")

outlier_map                = load_carrier_map("outlier_carrier_map.txt")
nonoutlier_map             = load_carrier_map("nonoutlier_carrier_map.txt")
parent_outlier_map         = load_carrier_map("parent_outlier_carrier_map.txt")
parent_nonoutlier_map      = load_carrier_map("parent_nonoutlier_carrier_map.txt")

fix_bed(
    "~{cohort_prefix}.outlier_samples.bed.gz",
    "~{cohort_prefix}.outlier_samples.fixed.bed.gz",
    outlier_map
)
fix_bed(
    "~{cohort_prefix}.nonoutlier_samples.bed.gz",
    "~{cohort_prefix}.nonoutlier_samples.fixed.bed.gz",
    nonoutlier_map
)
fix_bed(
    "~{cohort_prefix}.no_proband_carriers_outlier.bed.gz",
    "~{cohort_prefix}.no_proband_carriers_outlier.fixed.bed.gz",
    parent_outlier_map
)
fix_bed(
    "~{cohort_prefix}.no_proband_carriers_nonoutlier.bed.gz",
    "~{cohort_prefix}.no_proband_carriers_nonoutlier.fixed.bed.gz",
    parent_nonoutlier_map
)

print("Done.")
PYEOF

    # Replace originals with fixed versions
    mv "${PREFIX}.outlier_samples.fixed.bed.gz"                  "${PREFIX}.outlier_samples.bed.gz"
    mv "${PREFIX}.nonoutlier_samples.fixed.bed.gz"               "${PREFIX}.nonoutlier_samples.bed.gz"
    mv "${PREFIX}.no_proband_carriers_outlier.fixed.bed.gz"      "${PREFIX}.no_proband_carriers_outlier.bed.gz"
    mv "${PREFIX}.no_proband_carriers_nonoutlier.fixed.bed.gz"   "${PREFIX}.no_proband_carriers_nonoutlier.bed.gz"

    # Rename parent-only ID lists to final output names
    mv ids_parent_only_in_outliers.txt    ids_no_proband_carriers_outlier.txt
    mv ids_parent_only_in_nonoutliers.txt ids_no_proband_carriers_nonoutlier.txt

  >>>

  output {
    # Outlier-exclusive outputs (main PROBAND_AC branch)
    File outlier_vcf                              = "~{cohort_prefix}.outlier_samples.vcf.gz"
    File outlier_vcf_tbi                          = "~{cohort_prefix}.outlier_samples.vcf.gz.tbi"
    File outlier_bed                              = "~{cohort_prefix}.outlier_samples.bed.gz"

    # Non-outlier outputs (main PROBAND_AC branch, outlier samples fully excluded)
    File nonoutlier_vcf                           = "~{cohort_prefix}.nonoutlier_samples.vcf.gz"
    File nonoutlier_vcf_tbi                       = "~{cohort_prefix}.nonoutlier_samples.vcf.gz.tbi"
    File nonoutlier_bed                           = "~{cohort_prefix}.nonoutlier_samples.bed.gz"

    # All-samples BED (main PROBAND_AC branch)
    File all_samples_bed                          = "~{cohort_prefix}.all_samples.bed.gz"

    # Variant ID lists (main branch)
    File ids_outlier_exclusive                    = "ids_outlier_exclusive.txt"
    File ids_in_nonoutliers                       = "ids_in_nonoutliers.txt"

    # Parent-only outputs: PROBAND_AC=0 branch, split by outlier status
    File no_proband_carriers_outlier_vcf          = "~{cohort_prefix}.no_proband_carriers_outlier.vcf.gz"
    File no_proband_carriers_outlier_vcf_tbi      = "~{cohort_prefix}.no_proband_carriers_outlier.vcf.gz.tbi"
    File no_proband_carriers_outlier_bed          = "~{cohort_prefix}.no_proband_carriers_outlier.bed.gz"
    File no_proband_carriers_nonoutlier_vcf       = "~{cohort_prefix}.no_proband_carriers_nonoutlier.vcf.gz"
    File no_proband_carriers_nonoutlier_vcf_tbi   = "~{cohort_prefix}.no_proband_carriers_nonoutlier.vcf.gz.tbi"
    File no_proband_carriers_nonoutlier_bed       = "~{cohort_prefix}.no_proband_carriers_nonoutlier.bed.gz"
    File ids_no_proband_carriers_outlier          = "ids_no_proband_carriers_outlier.txt"
    File ids_no_proband_carriers_nonoutlier       = "ids_no_proband_carriers_nonoutlier.txt"
    File proband_ids                              = "proband_ids.txt"
  }
}
