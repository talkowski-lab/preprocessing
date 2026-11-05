#!/usr/bin/env python3
"""
Regression tests for scripts/hail_merge_vcf_samples_update_info.py, run on the small VCFs in data/.
Each test runs the merge script as the WDL does and checks the merged VCF.

    python3 tests/merge_vcf_samples/test_merge_vcf_samples.py [test_name ...]

Needs Hail. Takes a few minutes (one Hail run per merge).

data/ holds three "rounds" with overlapping alleles at chr1:1000 (min-rep after splitting):
  roundA  A1 A2 A3   AT>A,ATT   -> -T (AT>A) and +T (A>AT); A3 is 0/2
                     T>A at 2000, where A3 is a no-call (./.)
  roundB  B1 B2      AT>A,ATTT  -> -T (AT>A) and +TT (A>ATT); B2 is 0/2
                     G>C at 3000, FILTER LowQual, where B2 is a no-call (./.)
  roundC  C1 C2      A>AT       -> +T, written with a different REF than in roundA
                     G>C at 3000, FILTER PASS
  roundD_dup_A1      roundC's records with samples A1 (also in roundA) and D1
"""
import gzip
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
SCRIPT = os.path.join(REPO, 'scripts', 'hail_merge_vcf_samples_update_info.py')
DATA = os.path.join(HERE, 'data')

# variant keys after splitting: (pos, ref, alt)
DEL_T = (1000, 'AT', 'A')
INS_T = (1000, 'A', 'AT')
INS_TT = (1000, 'A', 'ATT')
SNV_2000 = (2000, 'T', 'A')
SNV_3000 = (3000, 'G', 'C')


def read_vcf(path):
    """Merged VCF as {(pos, ref, alt): {'qual', 'filter', 'info', 'format', 'gt': {sample: GT},
    'entries': {sample: {FORMAT field: value}}}} plus samples and FORMAT keys"""
    rows, samples, format_keys = {}, [], set()
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('##FORMAT=<ID='):
                format_keys.add(line.split('ID=')[1].split(',')[0])
            elif line.startswith('#CHROM'):
                samples = line.rstrip('\n').split('\t')[9:]
            elif not line.startswith('#'):
                fields = line.rstrip('\n').split('\t')
                info = dict(kv.split('=', 1) if '=' in kv else (kv, True) for kv in fields[7].split(';'))
                gts = {s: entry.split(':')[0] for s, entry in zip(samples, fields[9:])}
                # trailing missing FORMAT values may be left off an entry, so absent keys mean missing
                entries = {s: dict(zip(fields[8].split(':'), entry.split(':'))) for s, entry in zip(samples, fields[9:])}
                rows[(int(fields[1]), fields[3], fields[4])] = {
                    'qual': float(fields[5]), 'filter': fields[6], 'info': info,
                    'format': fields[8], 'gt': gts, 'entries': entries,
                }
    return rows, samples, format_keys


def run_merge(tmp, name, vcfs, *extra, expect_fail=False):
    out = os.path.join(tmp, f'{name}.vcf.bgz')
    cmd = [sys.executable, SCRIPT,
           '--vcf-files', *[os.path.join(DATA, v) for v in vcfs],
           '--batch-size', '2',  # 3+ inputs also exercise the recursive union
           '--output-vcf-file', out,
           '--mem', '4',
           '--min-info-fields', 'VQSLOD',
           '--sum-info-fields', 'DP',
           *extra]
    res = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True)
    log = res.stdout + res.stderr
    if expect_fail:
        assert res.returncode != 0, f"{name}: expected the merge to fail"
        return log
    assert res.returncode == 0, f"{name}: merge failed\n{log[-3000:]}"
    assert os.path.exists(out + '.tbi'), f"{name}: no tabix index"
    return read_vcf(out) + (log,)


def check_row(rows, key, ac=None, an=None, gt=None, **fields):
    assert key in rows, f"missing row {key}; rows: {sorted(rows)}"
    row = rows[key]
    if ac is not None:
        assert int(row['info']['AC']) == ac == int(row['info']['cohort_AC']), f"{key} AC {row['info']['AC']} / cohort_AC {row['info']['cohort_AC']}, expected {ac}"
    if an is not None:
        assert int(row['info']['AN']) == an == int(row['info']['cohort_AN']), f"{key} AN {row['info']['AN']} / cohort_AN {row['info']['cohort_AN']}, expected {an}"
    for sample, expected in (gt or {}).items():
        assert row['gt'][sample] == expected, f"{key} {sample}: GT {row['gt'][sample]}, expected {expected}"
    for field, expected in fields.items():
        actual = row[field] if field in ('qual', 'filter') else row['info'][field]
        assert str(actual) == str(expected), f"{key} {field}: {actual}, expected {expected}"


def comparable(rows):
    """Everything that should not depend on input order: GTs by sample, counts, QUAL, FILTER, min/sum fields.
    INFO fields merged by 'first non-missing' (e.g. AS_VQSLOD) take the first input's value, so they do."""
    return {key: (row['qual'], row['filter'], row['gt'],
                  {k: row['info'][k] for k in ('AC', 'AN', 'AF', 'VQSLOD', 'DP')})
            for key, row in rows.items()}


def test_default_merge(tmp):
    """outer join, fill_missing=absent: one row per allele, per-allele INFO, counts over all 7 samples"""
    rows, samples, _, _ = run_merge(tmp, 'default', ['roundA.vcf', 'roundB.vcf', 'roundC.vcf'])
    assert sorted(rows) == sorted([DEL_T, INS_T, INS_TT, SNV_2000, SNV_3000]), sorted(rows)
    # -T from roundA and roundB in one row; roundC has no -T, so its samples are filled 0/0
    check_row(rows, DEL_T, ac=4, an=14, gt={'A1': '0/1', 'A2': '1/1', 'A3': '0/0', 'B1': '0/1', 'B2': '0/0', 'C1': '0/0', 'C2': '0/0'},
              qual=40.0, VQSLOD='1.50000', DP=50, AS_VQSLOD='1.50000')
    # +T written as AT>ATT (roundA) and A>AT (roundC) in one row
    check_row(rows, INS_T, ac=2, an=14, gt={'A3': '0/1', 'C1': '0/1', 'A1': '0/0', 'B1': '0/0'}, qual=45.0, DP=55)
    check_row(rows, INS_TT, ac=1, an=14, gt={'B2': '0/1', 'A3': '0/0', 'C1': '0/0'})
    # genuine no-calls stay ./. (A3 at 2000, B2 at 3000); other rounds' absent samples are filled
    check_row(rows, SNV_2000, ac=1, an=12, gt={'A3': './.', 'B1': '0/0', 'C2': '0/0'})
    check_row(rows, SNV_3000, ac=2, an=12, gt={'B2': './.', 'A1': '0/0', 'C2': '0/1'}, filter='LowQual', qual=40.0, DP=45)


def test_input_order(tmp):
    """the same inputs in reverse order give the same genotypes, counts, QUAL, FILTER and min/sum fields"""
    forward = run_merge(tmp, 'order_abc', ['roundA.vcf', 'roundB.vcf', 'roundC.vcf'])[0]
    reverse = run_merge(tmp, 'order_cba', ['roundC.vcf', 'roundB.vcf', 'roundA.vcf'])[0]
    assert comparable(forward) == comparable(reverse), "merge depends on input order"


def test_fill_all(tmp):
    """fill_missing=all also turns genuine no-calls into 0/0"""
    rows = run_merge(tmp, 'fill_all', ['roundA.vcf', 'roundB.vcf', 'roundC.vcf'], '--fill-missing', 'all')[0]
    check_row(rows, SNV_2000, ac=1, an=14, gt={'A3': '0/0'})
    check_row(rows, SNV_3000, ac=2, an=14, gt={'B2': '0/0'})


def test_fill_none(tmp):
    """fill_missing=none: samples whose round lacks the variant stay missing"""
    rows = run_merge(tmp, 'fill_none', ['roundA.vcf', 'roundB.vcf', 'roundC.vcf'], '--fill-missing', 'none')[0]
    check_row(rows, DEL_T, ac=4, an=10, gt={'C1': './.'})
    check_row(rows, INS_T, ac=2, an=10, gt={'B1': './.'})
    check_row(rows, INS_TT, ac=1, an=4, gt={'A1': './.', 'C1': './.'})
    check_row(rows, SNV_2000, ac=1, an=4, gt={'A3': './.', 'B1': './.'})
    check_row(rows, SNV_3000, ac=2, an=6, gt={'A1': './.'})


def test_fill_absent_at_genotyped_loci(tmp):
    """fill only where the sample's round has a record at the same locus (every round has chr1:1000)"""
    rows = run_merge(tmp, 'fill_loci', ['roundA.vcf', 'roundB.vcf', 'roundC.vcf'], '--fill-missing', 'absent_at_genotyped_loci')[0]
    check_row(rows, DEL_T, ac=4, an=14, gt={'C1': '0/0'})
    check_row(rows, INS_TT, ac=1, an=14, gt={'A1': '0/0', 'C1': '0/0'})
    check_row(rows, SNV_2000, ac=1, an=4, gt={'A3': './.', 'B1': './.', 'C1': './.'})
    check_row(rows, SNV_3000, ac=2, an=6, gt={'A1': './.', 'B2': './.', 'C1': '0/0'})


def test_rows_in_one_input(tmp):
    """a row present in only one input keeps that input's QUAL, FILTER and INFO, whichever side of
    union_cols it comes from (other inputs' copies of these fields are missing for that row)"""
    # one batch: roundA is the left MT, roundB the right
    rows = run_merge(tmp, 'one_input_ab', ['roundA.vcf', 'roundB.vcf'])[0]
    check_row(rows, INS_T, qual=50.0, filter='PASS', VQSLOD='1.50000', DP=30, AS_VQSLOD='-2.50000')  # roundA only
    check_row(rows, SNV_2000, qual=50.0, filter='PASS', VQSLOD='3.00000', DP=30, AS_VQSLOD='3.00000')  # roundA only
    check_row(rows, INS_TT, qual=40.0, filter='PASS', VQSLOD='2.00000', DP=20, AS_VQSLOD='1.00000')  # roundB only
    check_row(rows, SNV_3000, qual=40.0, filter='LowQual', VQSLOD='4.00000', DP=20, AS_VQSLOD='4.00000')  # roundB only
    # batch size 2: roundB is alone in the second batch, so its rows come in through the pairwise union
    rows = run_merge(tmp, 'one_input_acb', ['roundA.vcf', 'roundC.vcf', 'roundB.vcf'])[0]
    check_row(rows, INS_TT, ac=1, an=14, qual=40.0, filter='PASS', VQSLOD='2.00000', DP=20, AS_VQSLOD='1.00000')
    check_row(rows, SNV_2000, qual=50.0, filter='PASS', VQSLOD='3.00000', DP=30)


def test_haploid_calls_diploidized(tmp):
    """haploid calls (DRAGEN males, non-PAR chrX) become homozygous diploid before splitting;
    their PL is dropped; diploid calls and their PL are unchanged; AN counts every sample as diploid"""
    rows = run_merge(tmp, 'haploid', ['haploidA_chrX.vcf', 'diploidB_chrX.vcf'])[0]
    snv_a, snv_t, snv_20m = (10000000, 'G', 'A'), (10000000, 'G', 'T'), (20000000, 'C', 'T')
    assert sorted(rows) == sorted([snv_a, snv_t, snv_20m]), sorted(rows)
    # H1 is haploid 2 at G>A,T -> 2/2 -> 0/0 for G>A and 1/1 for G>T
    check_row(rows, snv_a, ac=1, an=6, gt={'H1': '0/0', 'H2': '0/1', 'K1': '0/0'})
    check_row(rows, snv_t, ac=2, an=6, gt={'H1': '1/1', 'H2': '0/0', 'K1': '0/0'})
    # H1 is haploid 0 at 20M -> 0/0
    check_row(rows, snv_20m, ac=2, an=6, gt={'H1': '0/0', 'H2': '0/1', 'K1': '0/1'})
    for key in rows:
        assert rows[key]['entries']['H1'].get('PL', '.') == '.', f"{key} H1 PL should be missing: {rows[key]['entries']['H1']}"
    assert rows[snv_20m]['entries']['H2']['PL'] == '50,0,50', rows[snv_20m]['entries']['H2']
    assert rows[snv_a]['entries']['H1']['AD'] == '10,0', rows[snv_a]['entries']['H1']  # AD still downcoded


def test_inner_join(tmp):
    """site_join=inner keeps only variants in every VCF"""
    rows = run_merge(tmp, 'inner', ['roundA.vcf', 'roundB.vcf'], '--site-join', 'inner')[0]
    assert sorted(rows) == [DEL_T], sorted(rows)
    check_row(rows, DEL_T, ac=4, an=10)


def test_subset_rename_drop(tmp):
    """sample file + cohort prefixes + format_fields_to_drop; a VCF left with no samples is skipped"""
    rows, samples, format_keys, log = run_merge(
        tmp, 'subset', ['roundA.vcf', 'roundB.vcf', 'roundC.vcf', 'roundD_dup_A1.vcf'],
        '--sample-file', os.path.join(DATA, 'samples_A1_B2_X9.txt'),
        '--cohort-prefixes', 'rA', 'rB', 'rC', 'rD',
        '--format-fields-to-drop', 'PL')
    assert samples == ['rA:A1', 'rB:B2', 'rD:A1'], samples
    assert 'Skipping' in log and 'roundC.vcf' in log, "roundC (no kept samples) should be skipped"
    assert "X9" in log, "sample X9 (in no VCF) should be reported"
    assert 'PL' not in format_keys and {'GT', 'AD'} <= format_keys, format_keys
    assert all(row['format'] == 'GT:AD' or row['format'] == 'GT' for row in rows.values())
    check_row(rows, DEL_T, ac=1, an=6, gt={'rA:A1': '0/1', 'rB:B2': '0/0', 'rD:A1': '0/0'})
    check_row(rows, INS_T, ac=1, an=6, gt={'rD:A1': '0/1', 'rA:A1': '0/0'})


def test_duplicate_samples_fail(tmp):
    """the same sample ID in two VCFs stops the merge unless prefixes make the IDs unique"""
    log = run_merge(tmp, 'dup', ['roundA.vcf', 'roundD_dup_A1.vcf'], expect_fail=True)
    assert 'appear in more than one VCF' in log, log[-2000:]


def test_drop_gt_fails(tmp):
    log = run_merge(tmp, 'drop_gt', ['roundA.vcf'], '--format-fields-to-drop', 'GT', expect_fail=True)
    assert 'GT cannot be dropped' in log, log[-2000:]


def test_single_input(tmp):
    """one VCF (the WDL's merge_inputs=false): split, per-allele INFO, counts recomputed, no-calls kept"""
    rows = run_merge(tmp, 'single', ['roundA.vcf'])[0]
    assert sorted(rows) == sorted([DEL_T, INS_T, SNV_2000]), sorted(rows)
    check_row(rows, DEL_T, ac=3, an=6, AS_VQSLOD='1.50000')
    check_row(rows, INS_T, ac=1, an=6, AS_VQSLOD='-2.50000')
    check_row(rows, SNV_2000, ac=1, an=4, gt={'A3': './.'})


TESTS = [test_default_merge, test_input_order, test_fill_all, test_fill_none,
         test_fill_absent_at_genotyped_loci, test_inner_join, test_subset_rename_drop,
         test_duplicate_samples_fail, test_drop_gt_fails, test_single_input, test_rows_in_one_input,
         test_haploid_calls_diploidized]


def main():
    selected = [t for t in TESTS if not sys.argv[1:] or t.__name__ in sys.argv[1:]]
    failed = []
    for test in selected:
        with tempfile.TemporaryDirectory() as tmp:
            try:
                test(tmp)
                print(f"PASS {test.__name__}", flush=True)
            except AssertionError as e:
                failed.append(test.__name__)
                print(f"FAIL {test.__name__}: {e}", flush=True)
    print(f"\n{len(selected) - len(failed)}/{len(selected)} passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
