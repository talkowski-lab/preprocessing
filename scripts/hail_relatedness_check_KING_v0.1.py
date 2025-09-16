###
# DESCRIPTION: TODO

## CHANGE LOG:
'''
1/31/2025:
- added split_multi parameter
- removed somalier_vcf input (same as sites_uri in the WDL)
'''
###

import pandas as pd
import numpy as np
import hail as hl
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import ast
import datetime
from gnomad.resources.grch38 import gnomad
from gnomad.sample_qc import relatedness

vcf_uri = sys.argv[1]
cohort_prefix = sys.argv[2]
ped_uri = sys.argv[3]
cores = sys.argv[4]  # string
mem = int(np.floor(float(sys.argv[5])))
bucket_id = sys.argv[6]
score_table = sys.argv[7]  # TODO: REMOVE THIS OUTPUT (only here for compatibility with PC-Relate WDL task)
genome_build = sys.argv[8]
split_multi = ast.literal_eval(sys.argv[9].capitalize())
kinship_ht_uri = sys.argv[10]  # optional

hl.init(min_block_size=128, 
        local=f"local[*]", 
        spark_conf={
                    "spark.driver.memory": f"{int(np.floor(mem*0.8))}g",
                    "spark.speculation": 'true'
                    }, 
        tmp_dir="tmp", local_tmpdir="tmp",
                    )

#split-multi
def split_multi_ssc(mt):
    mt = mt.annotate_rows(num_alleles = mt.alleles.size() ) # Add number of alleles at site before split
    # only split variants that aren't already split
    bi = mt.filter_rows(hl.len(mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=1, was_split=False, old_locus=bi.locus, old_alleles=bi.alleles)
    multi = mt.filter_rows(hl.len(mt.alleles) > 2)
    # Now split
    split = hl.split_multi(multi, permit_shuffle=True)
    sm = split.union_rows(bi)
    # sm = hl.split_multi(mt, permit_shuffle=True)
    if 'PL' in list(mt.entry.keys()):
        pl = hl.or_missing(hl.is_defined(sm.PL),
                        (hl.range(0, 3).map(lambda i: hl.min(hl.range(0, hl.len(sm.PL))
        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j), sm.a_index) == hl.unphased_diploid_gt_index_call(i))
        .map(lambda j: sm.PL[j])))))
        sm = sm.annotate_entries(PL = pl)
    split_ds = sm.annotate_entries(GT = hl.downcode(sm.GT, sm.a_index),
                                   AD = hl.or_missing(hl.is_defined(sm.AD), [hl.sum(sm.AD) - sm.AD[sm.a_index], sm.AD[sm.a_index]])
                                   ) 
        #GQ = hl.cond(hl.is_defined(pl[0]) & hl.is_defined(pl[1]) & hl.is_defined(pl[2]), hl.gq_from_pl(pl), sm.GQ) )
    mt = split_ds.drop('old_locus', 'old_alleles')
    return mt

mt = hl.import_vcf(vcf_uri, reference_genome=genome_build, force_bgz=True, call_fields=[], array_elements_required=False)
if split_multi:
    mt = split_multi_ssc(mt)

# RUN KING AND IDENTITY_BY_DESCENT
king_mt = hl.king(mt.GT)
king_ht = king_mt.entries().rename({'s_1':'i', 's':'j'}).key_by('i','j')

ibd_ht = hl.identity_by_descent(mt)

rel = ibd_ht.annotate(phi=king_ht[ibd_ht.key].phi).flatten().key_by('i','j')

# ANNOTATE RELATIONSHIP FIELD USING KING HARD CUTOFFS
# SOURCE 1: https://www.cog-genomics.org/plink/2.0/distance#king_coefs
# Simple KING cutoffs for inferred relatedness
# Note that KING kinship coefficients are scaled such that duplicate samples have kinship 0.5, not 1. 
# First-degree relations (parent-child, full siblings) correspond to ~0.25, 
# second-degree relations correspond to ~0.125, etc. 
# It is conventional to use a cutoff of ~0.354 (the geometric mean of 0.5 and 0.25) 
# to screen for monozygotic twins and duplicate samples, ~0.177 to add first-degree relations, etc.

# SOURCE 2: https://bioweb.pasteur.fr/docs/modules/king/1.4/#:~:text=Parameter%20%2D%2Drelated%20%2D%2Ddegree,be%20excluded%20from%20the%20output.
# lose relatives can be inferred fairly reliably based on the estimated kinship coefficients 
# as shown in the following simple algorithm: an estimated kinship coefficient range 
# >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] corresponds to 
# duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively. 

# can't distinguish between parent-child and siblings just using kinship coefficient (fine for now?)
# king_cutoffs = {'parent-child': [0.177, 0.354],  
#                'second degree relatives': [0.0884, 0.177],
#                'duplicate/twins': [0.354, 1],
#                'unrelated': [0, 0.0884]}

# rel = rel.annotate(relationship=hl.case()
#                  .when(rel.phi >= king_cutoffs['duplicate/twins'][0], 'duplicate/twins')
#                  .when(rel.phi >= king_cutoffs['parent-child'][0], 'parent-child')
#                  .when(rel.phi >= king_cutoffs['second degree relatives'][0], 'second degree relatives')
#                  .when(rel.phi >= king_cutoffs['unrelated'][0], 'unrelated')
#                  .or_missing())

# ANNOTATE RELATIONSHIP USING PC-RELATE CUTOFFS
rel = rel.annotate(relationship=relatedness.get_relationship_expr(kin_expr = rel.phi, ibd0_expr = rel['ibd.Z0'], 
                                                            ibd1_expr = rel['ibd.Z1'], ibd2_expr = rel['ibd.Z2']))

# Optionally write HT
if kinship_ht_uri!='NA':
    rel = rel.checkpoint(kinship_ht_uri, overwrite=True)

ped = pd.read_csv(ped_uri, sep='\t').iloc[:, :6]
ped.columns = ['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'phenotype']
ped.index = ped.sample_id
try:
    ped['sex'] = ped.sex.astype(str).str.lower().replace({'male': 1, 'female': 2})
    ped.loc[~ped.sex.isin([1,2,'1','2']), 'sex'] = 0
    ped['sex'] = ped.sex.astype(int)
except Exception as e:
    print(str(e))
    pass
ped[['family_id', 'sample_id', 'paternal_id', 'maternal_id']] = ped[['family_id', 'sample_id', 'paternal_id', 'maternal_id']].astype(str)

# subset ped to samples in VCF
vcf_samps = mt.s.collect()
ped = ped[ped.sample_id.isin(vcf_samps)].copy()

fam_sizes = ped.family_id.value_counts().to_dict()
fathers = ped[~ped.paternal_id.isin(['0','-9'])].set_index('sample_id').paternal_id.to_dict()
mothers = ped[~ped.maternal_id.isin(['0','-9'])].set_index('sample_id').maternal_id.to_dict()

def get_sample_role(row):
    if fam_sizes[row.family_id]==1:
        role = 'Singleton'
    elif (row.maternal_id=='0') & (row.paternal_id=='0'):
        if (row.sample_id in fathers.values()):
            role = 'Father'
        elif (row.sample_id in mothers.values()):
            role = 'Mother'
        else:
            role = 'Unknown'
    elif (row.maternal_id in mothers.values()) & (row.paternal_id in fathers.values()):
        role = 'Proband'
    else:
        role = 'Unknown'
    return role

if not ped.empty:
    ped['role'] = ped.apply(get_sample_role, axis=1)

    ped_ht = hl.Table.from_pandas(ped)

    dad_temp = ped_ht.key_by('sample_id','paternal_id')
    dad_temp2 = ped_ht.key_by('paternal_id','sample_id')

    all_dad_ht = dad_temp.join(rel.semi_join(dad_temp)).key_by().union(dad_temp2.join(rel.semi_join(dad_temp2)).key_by(), unify=True)
    all_dad_ht = all_dad_ht.annotate(father_status=all_dad_ht.relationship).drop('relationship')

    mom_temp = ped_ht.key_by('sample_id','maternal_id')
    mom_temp2 = ped_ht.key_by('maternal_id','sample_id')

    all_mom_ht = mom_temp.join(rel.semi_join(mom_temp)).key_by().union(mom_temp2.join(rel.semi_join(mom_temp2)).key_by(), unify=True)
    all_mom_ht = all_mom_ht.annotate(mother_status=all_mom_ht.relationship).drop('relationship')

    mom_df = all_mom_ht.to_pandas()
    dad_df = all_dad_ht.to_pandas()

    # CHANGED FOR KING AND IDENTITY_BY_DESCENT
    rename_cols = list(np.setdiff1d(list(rel.row), ['i', 'j']))
    mom_df = mom_df.rename({col: 'mother_'+col for col in rename_cols}, axis=1).copy()
    dad_df = dad_df.rename({col: 'father_'+col for col in rename_cols}, axis=1).copy()

    all_df = mom_df.merge(dad_df, how='outer')

    merged_rel_df = pd.concat([ped, all_df.set_index('sample_id')], axis=1)
    merged_rel_df = merged_rel_df.loc[:,~merged_rel_df.columns.duplicated()].copy()

else:
    merged_rel_df = ped.copy()

merged_rel_df.to_csv(f"{cohort_prefix}_relatedness_qc.ped", sep='\t', index=False)

# annotate ped
ped_rels = {'i':[], 'j': [], 'ped_relationship': [], 'family_id': []}

for fam in ped.family_id.unique():
    fam_df = ped[ped.family_id==fam].reset_index(drop=True)
    for i, row_i in fam_df.iterrows():
        for j, row_j in fam_df.iterrows():
            if (i==j) | ((row_i.role in ['Mother','Father']) & (row_j.role in ['Mother','Father'])):
                continue

            if (((row_i.paternal_id == row_j.sample_id)) | ((row_i.sample_id == row_j.paternal_id))) |\
            (((row_i.maternal_id == row_j.sample_id)) | ((row_i.sample_id == row_j.maternal_id))):                
                ped_rels['ped_relationship'].append('parent-child')
            elif (row_i.paternal_id==row_j.paternal_id) & (row_i.maternal_id==row_j.maternal_id):
                    ped_rels['ped_relationship'].append('siblings')
            else:
                ped_rels['ped_relationship'].append('related_other')      

            ped_rels['i'].append(row_i.sample_id)
            ped_rels['j'].append(row_j.sample_id)
            ped_rels['family_id'].append(fam)


ped_rels_df = pd.DataFrame(ped_rels)
ped_rels_ht = hl.Table.from_pandas(ped_rels_df)

ped_rels_ht_merged = ped_rels_ht.annotate(i=ped_rels_ht.j, j=ped_rels_ht.i).key_by('i', 'j').union(ped_rels_ht.key_by('i','j'))

rel_merged = rel.key_by()
rel_merged = rel_merged.annotate(i=rel_merged.j, j=rel_merged.i).key_by('i', 'j').union(rel.key_by('i','j'))

try:
    related_in_ped = rel_merged.annotate(ped_relationship=ped_rels_ht_merged[rel_merged.key].ped_relationship)
except:  # no related in ped
    related_in_ped = rel_merged.annotate(ped_relationship=hl.missing('str'))
related_in_ped = related_in_ped.filter(hl.is_missing(related_in_ped.ped_relationship), keep=False)

unrelated_in_ped = rel_merged.anti_join(related_in_ped).annotate(ped_relationship='unrelated')

p = 0.05
only_related = unrelated_in_ped.filter(unrelated_in_ped.relationship!='unrelated')
downsampled_unrelated = unrelated_in_ped.filter(unrelated_in_ped.relationship=='unrelated').sample(p)

rel_total = related_in_ped.union(only_related).union(downsampled_unrelated)

rel_df = rel_total.to_pandas()
rel_df.to_csv(f"{cohort_prefix}_kinship.tsv.gz", sep='\t', index=False)