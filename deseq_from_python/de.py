#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# # Gene differential expression with DESeq2 from Python
#
# Author: Charles E. Vejnar (charles.vejnar@gmail.com)
#
# **Requirements**
#
# * Python with [Matplotlib](https://matplotlib.org)
# * Download gene expression data [dev_timecourse_zebrafish_vejnar_giraldez_count_900_gene_protein_coding_v11.csv.xz](https://data.giraldezlab.org/pub/vejnar_et_al_2019_genome_research/dev_timecourse_gene_count/dev_timecourse_zebrafish_vejnar_giraldez_count_900_gene_protein_coding_v11.csv.xz). More samples are [available](https://www.giraldezlab.org/data/vejnar_et_al_2019_genome_research/).
#
# Please cite relevant papers if you use this tutorial and/or data:
# * [Genome wide analysis of 3' UTR sequence elements and proteins regulating mRNA stability during maternal-to-zygotic transition in zebrafish](https://pubmed.ncbi.nlm.nih.gov/31227602/) Vejnar et al, 2019
# * [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2](https://pubmed.ncbi.nlm.nih.gov/25516281/) Love et al, 2014
#

import subprocess

import numpy as np
import pandas as pd

# Path to gene count
path_count = 'dev_timecourse_zebrafish_vejnar_giraldez_count_900_gene_protein_coding_v11.csv.xz'
# Minimum read count per condition
min_count = 1
# Method for adjusting P-values
p_adjust = 'fdr'
# Column separator
col_sep = ' '
# Define DE tests
tests = [{'name': 'WT shield vs a-Am',
          'samples1': ['WT shield pA B1', 'WT shield pA B2'],
          'samples2': ['WT a-Am shield pA B1', 'WT a-Am shield pA B2']}]

# Open main time-course
dcount = pd.read_csv(path_count, index_col=0)
# Keep only zebrafish
dcount.drop([i for i in dcount.index if not i.startswith('ENSDAR')], inplace=True)

de_data = [dcount.loc[:,['gene_name', 'gene_length']]]
for test in tests:
    # Get samples & Remove total
    dcount_cond = dcount.loc[:, test['samples1']+test['samples2']]
    # Selection min_count per condition
    sel_count = np.all([np.any(dcount_cond.loc[:, test['samples1']] >= min_count, axis=1), np.any(dcount_cond.loc[:, test['samples2']] >= min_count, axis=1)], axis=0)
    # Round (DESeq2 requires counts i.e. integers) & Output
    dcount_cond.loc[sel_count, :].round().to_csv('count.csv', index=True)

    # Conditions
    pd.DataFrame([test['samples1']+test['samples2'], ['a']*len(test['samples1']) + ['b']*len(test['samples2'])], index=['', 'condition']).T.to_csv('cond.csv', index=False)

    # Run DESeq
    subprocess.run(['Rscript','run_deseq.r', 'count.csv', 'cond.csv', 'de.csv', p_adjust], check=True)

    # Get results
    de = pd.read_csv('de.csv', index_col=0)
    de.rename(columns=lambda x: test['name']+col_sep+x, inplace=True)
    de_data.append(de)

# Merge and save
dds = pd.concat(de_data, join='inner', axis=1)
dds.to_csv('de.csv')
