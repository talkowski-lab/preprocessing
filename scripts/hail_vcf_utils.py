"""
Shared Hail helpers for VCF processing.

Fetched next to the calling script at runtime (see the hail_vcf_utils_script WDL inputs), so the
caller can `from hail_vcf_utils import ...`. Pin the URL to a tag or commit, not a branch:
every pipeline that splits alleles must split them the same way for joins on (locus, alleles).
"""
import hail as hl


def split_multi_ssc(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Split multiallelic sites. Biallelic sites are left as they are; split alleles are minimally
    represented by hl.split_multi. Adds row fields num_alleles, a_index and was_split.
    GT, AD and PL are downcoded to the split allele. INFO is not changed: call
    subset_allele_info_fields() afterwards.
    """
    mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
    # only split variants that aren't already split
    bi = mt.filter_rows(hl.len(mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=1, was_split=False, old_locus=bi.locus, old_alleles=bi.alleles)
    multi = mt.filter_rows(hl.len(mt.alleles) > 2)
    # Now split
    split = hl.split_multi(multi, permit_shuffle=True)
    sm = split.union_rows(bi)
    if 'PL' in list(mt.entry.keys()):
        pl = hl.or_missing(hl.is_defined(sm.PL),
                        (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
        .map(lambda j: sm.PL[j])))))
        sm = sm.annotate_entries(PL = pl)
    if 'AD' in list(mt.entry.keys()):
        sm = sm.annotate_entries(AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]]))
    split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index))
    mt = split_ds.drop('old_locus', 'old_alleles')
    return mt


def subset_allele_info_fields(mt: hl.MatrixTable, info_header: dict) -> hl.MatrixTable:
    """
    Subset per-allele INFO arrays (Number=A/R in info_header, i.e. hl.get_vcf_metadata(vcf)['info'])
    to the split allele, so split rows match their header. Requires a_index from split_multi_ssc().
    """
    for field, field_header in info_header.items():
        if (field not in list(mt.info)) or (not isinstance(mt.info[field].dtype, hl.tarray)):
            continue
        values = mt.info[field]
        if field_header['Number'] == 'A':
            mt = mt.annotate_rows(info=mt.info.annotate(**{field: hl.or_missing(hl.len(values) >= mt.a_index,
                                                                                [values[mt.a_index - 1]])}))
        elif field_header['Number'] == 'R':
            mt = mt.annotate_rows(info=mt.info.annotate(**{field: hl.or_missing(hl.len(values) > mt.a_index,
                                                                                [values[0], values[mt.a_index]])}))
    return mt
