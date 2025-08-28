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

import argparse
import os

parser = argparse.ArgumentParser(description="Convert BED-like input to VCF with defined row/col fields.")

parser.add_argument("--bed-uri", required=True, help="Input BED file URI")
parser.add_argument("--row-key", default="", help="Comma-separated list of row key fields")
parser.add_argument("--col-key", default="", help="Comma-separated list of column key fields")
parser.add_argument("--skip-fields", default="", help="Comma-separated list of fields to skip")
parser.add_argument("--col-fields", default="", help="Comma-separated list of column fields")
parser.add_argument("--entry-fields", default="", help="Comma-separated list of entry fields")
parser.add_argument("--priority-row-fields", default="", help="Comma-separated list of priority row fields")
parser.add_argument("--genome-build", default="hg38", help="Genome build (e.g. hg19, hg38)")

args = parser.parse_args()

bed_uri = args.bed_uri
row_key = args.row_key.split(",") if args.row_key else []
col_key = args.col_key.split(",") if args.col_key else []
skip_fields = (args.skip_fields.split(",") if args.skip_fields else []) + row_key
col_fields = args.col_fields.split(",") if args.col_fields else []
entry_fields = args.entry_fields.split(",") if args.entry_fields else []
priority_row_fields = args.priority_row_fields.split(",") if args.priority_row_fields else []
genome_build = args.genome_build

file_ext = bed_uri.split(".")[-1]
output_filename = os.path.basename(bed_uri).split(".bed")[0] + ".vcf.bgz"

hl.init(default_reference=genome_build)

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
bed_ht = bed_ht.annotate(rsid=bed_ht.variant_name)

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
og_row_fields = list(np.setdiff1d(list(bed_mt.row), row_key + ['locus','alleles']))
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

# KEY ROWS BY LOCUS AND ALLELES FOR EXPORT
# TAKES A LOT LONGER WITH THE RSID FIXES ABOVE?
hl.export_vcf(bed_mt.key_rows_by('locus', 'alleles'), output_filename)