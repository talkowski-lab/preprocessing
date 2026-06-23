version 1.0

workflow count_hq_het_sites {

    input {
        File   vcf
        File   vcf_index
        Int    min_gq       = 20
        Int    min_ad_total = 10
        Float  min_ab       = 0.2
        Float  max_ab       = 0.8
        String docker
        Int    mem_gb       = 16
        Int    cpu          = 4
    }

    call CountHQHet {
        input:
            vcf          = vcf,
            vcf_index    = vcf_index,
            min_gq       = min_gq,
            min_ad_total = min_ad_total,
            min_ab       = min_ab,
            max_ab       = max_ab,
            docker       = docker,
            mem_gb       = mem_gb,
            cpu          = cpu
    }

    output {
        File hq_het_counts = CountHQHet.hq_het_counts
    }

    meta {
        author      : "Generated WDL"
        description : "Count HQ heterozygous SNV/Indel sites per chromosome per sample"
    }
}

task CountHQHet {
    input {
        File   vcf
        File   vcf_index
        Int    min_gq
        Int    min_ad_total
        Float  min_ab
        Float  max_ab
        String docker
        Int    mem_gb
        Int    cpu
        Int    disk_gb = ceil(size(vcf, "GiB") * 2) + 500
    }

    command <<<
        set -euo pipefail

        # Debug: show what FILTER values actually look like to bcftools
        echo "Sample FILTER values seen by bcftools:" >&2
        bcftools query -f '%FILTER\n' ~{vcf} | sort | uniq -c | sort -rn | head -10 >&2

        # Keep only PASS sites (handles both VQSR-filtered VCFs with "PASS"
        # and unfiltered GATK VCFs with ".")
        bcftools view \
            --exclude-type other \
            --include 'FILTER="PASS" || FILTER="."' \
            --threads ~{cpu} \
            ~{vcf} \
        | bcftools query \
            -f '[%CHROM\t%SAMPLE\t%GT\t%GQ\t%AD\n]' \
            -o raw_genotypes.tsv

        echo "Total lines extracted: $(wc -l < raw_genotypes.tsv)" >&2
        echo "First 5 lines:" >&2
        head -5 raw_genotypes.tsv >&2

        python3 -c "
import re, sys
from collections import defaultdict

min_gq       = ~{min_gq}
min_ad_total = ~{min_ad_total}
min_ab       = ~{min_ab}
max_ab       = ~{max_ab}

HET_RE = re.compile(r'^(0[/|][1-9][0-9]*|[1-9][0-9]*[/|]0)$')

counts = defaultdict(int)
skipped_gt = skipped_gq = skipped_ad = passed = 0

with open('raw_genotypes.tsv') as fh:
    for line in fh:
        parts = line.rstrip('\n').split('\t')
        if len(parts) < 5:
            continue
        chrom, sample, gt, gq_s, ad_s = parts

        if not HET_RE.match(gt):
            skipped_gt += 1
            continue

        try:
            if float(gq_s) < min_gq:
                skipped_gq += 1
                continue
        except ValueError:
            skipped_gq += 1
            continue

        try:
            ad     = ad_s.split(',')
            ref_ad = float(ad[0])
            alt_ad = float(ad[1])
            total  = ref_ad + alt_ad
            if total < min_ad_total:
                skipped_ad += 1
                continue
            if not (min_ab <= alt_ad / total <= max_ab):
                skipped_ad += 1
                continue
        except (ValueError, IndexError):
            skipped_ad += 1
            continue

        counts[(sample, chrom)] += 1
        passed += 1

print(f'Passed={passed} Skipped_GT={skipped_gt} Skipped_GQ={skipped_gq} Skipped_AD={skipped_ad}', file=sys.stderr)

with open('hq_het_counts.tsv', 'w') as out:
    out.write('sample\tchrom\tn_hq_het\n')
    for (sample, chrom), n in sorted(counts.items()):
        out.write(f'{sample}\t{chrom}\t{n}\n')

print(f'Done. Wrote {len(counts)} sample-chrom combinations.', file=sys.stderr)
"

    >>>

    output {
        File hq_het_counts = "hq_het_counts.tsv"
    }

    runtime {
        docker      : docker
        memory      : "~{mem_gb} GiB"
        cpu         : cpu
        disks       : "local-disk ~{disk_gb} HDD"
        preemptible : 2
    }
}
