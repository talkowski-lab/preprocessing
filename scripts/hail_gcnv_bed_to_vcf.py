###
# DESCRIPTION: TODO

## CHANGE LOG:
'''
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

bed_uri = sys.argv[1]
row_key = sys.argv[2].split(',')
col_key = sys.argv[3].split(',')
skip_fields = sys.argv[4].split(',') + row_key  # Also skip row_key when defining row_fields
col_fields = sys.argv[5].split(',')
entry_fields = sys.argv[6].split(',')
priority_row_fields = sys.argv[7].split(',')

file_ext = bed_uri.split('.')[-1]
output_filename = os.path.basename(bed_uri).split('.bed')[0] + '.vcf.bgz'

bed_df = pd.read_csv(bed_uri, sep='\t',compression='gzip' if file_ext in ['gz', 'bgz'] else None)

# Get all fields with semi-colons because they mess up VCF formatting
semi_colon_fields = []
for field in bed_df.columns:
    if bed_df[field].dtype!='object':
        continue
    if (bed_df[field].astype(str).str.contains(';').any()):
        semi_colon_fields.append(field)

field_types = {'start': 'int', 
                'end': 'int',
                'GT': 'int', 
                'CN': 'int',  
                'NP': 'int',  
                'QA': 'int',  
                'QS': 'int',  
                'QSE': 'int',  
                'QSS': 'int',  
                'ploidy': 'int', 
                'sc': 'int', 
                'sf': 'float'}

bed_ht = hl.import_table(bed_uri, force_bgz=file_ext=='bgz', force=file_ext!='bgz', types=field_types)  # cast start/end columns to int
bed_ht = bed_ht.annotate(start=bed_ht.start + 1)  # adjust for bed 0-based coordinates
bed_ht = bed_ht.annotate(alleles=['N', '<' + bed_ht.svtype + '>'],
                        locus=hl.locus(bed_ht.chr, bed_ht.start)).drop('chr','start')

# NEED UNIQUE ROW KEY TO NOT MESS THINGS UP
# Set variant_name field as rsid for VCF column with _REP{i} suffix for repeats
# Count occurrences of each variant_name
variant_name_count_dict = bed_ht.aggregate(hl.agg.counter(bed_ht.variant_name))
variant_name_count = hl.dict(variant_name_count_dict)

# Annotate rows with a new variant_name, adding a suffix for repeats
bed_ht = bed_ht.add_index()
bed_ht = bed_ht.annotate(
    rsid = hl.case()
        .when(variant_name_count[bed_ht.variant_name] > 1, 
              bed_ht.variant_name + "_REP" + hl.str(bed_ht.idx))
        .default(bed_ht.variant_name)
)

row_fields = list(np.setdiff1d(list(bed_ht.row), entry_fields + skip_fields + col_fields + col_key))

# Convert BED HT to MT for VCF export
bed_mt = bed_ht.to_matrix_table(row_key=row_key, 
                      col_key=col_key,
                      row_fields=row_fields,
                      col_fields=col_fields)

# Replace semi-colons
bed_mt = bed_mt.annotate_rows(**{field: bed_mt[field].replace(';',',') 
                        for field in semi_colon_fields})

# Rename end and svtype to match GATK-SV formatting (all caps)
bed_mt = bed_mt.annotate_rows(END=bed_mt.end, SVTYPE=bed_mt.svtype)
# Calculate SVLEN from start and end
bed_mt = bed_mt.annotate_rows(SVLEN=(bed_mt.end-bed_mt.locus.position)+1)
# Move row fields to INFO and drop original row fields
og_row_fields = list(np.setdiff1d(list(bed_mt.row), row_key))
bed_mt = bed_mt.annotate_rows(info=bed_mt.row.drop(*['end','svtype','locus','alleles'])).drop(*og_row_fields)

# Change order of rows
row_field_order = priority_row_fields + list(np.setdiff1d(list(bed_mt.info), priority_row_fields))
bed_mt = bed_mt.annotate_rows(info=bed_mt.info.select(*row_field_order))

# Count missing as hom_ref 
# Rename GT entry field from gCNV as gCNV_GT
bed_mt = bed_mt.annotate_entries(gCNV_GT=bed_mt.GT,
                                 GT=hl.if_else(hl.is_missing(bed_mt.GT), 
                                              hl.call(0, 0),
                                              hl.call(0, 1)))

variant_name_count_series = pd.Series(variant_name_count_dict)

# KEY ROWS BY LOCUS AND ALLELES FOR EXPORT
# TAKES A LOT LONGER WITH THE RSID FIXES ABOVE?
hl.export_vcf(bed_mt.key_rows_by('locus', 'alleles'), output_filename)