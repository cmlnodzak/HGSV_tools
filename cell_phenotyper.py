#!/usr/bin/py
import os
import pandas as pd
import numpy as np
import pybedtools
from pybedtools import BedTool


df = pd.read_table('/nobackup/shilab/Data/novogene/Counts/R697.counts')
df = pd.DataFrame(df)

gtf = BedTool('/nobackup/shilab/Data/novogene/Counts/gencode.v25.annotation.gtf').to_dataframe()

transcripts = gtf[gtf.feature == 'transcript']

annot = pd.DataFrame()
annot['gene_id'] = transcripts['attributes'].apply(lambda x: x.split('"')[1])
annot['gene_id'] = annot['gene_id'].apply(lambda x: x.split('.')[0])

annot['target_id'] = transcripts['attributes'].apply(lambda x: x.split('"')[3])
annot['target_id'] = annot['target_id'].apply(lambda x: x.split('.')[0])

annot = annot.reset_index().drop('index',1)

df_genes = df.merge(annot, on='target_id').drop('target_id', 1)

df_genes_group = df_genes.groupby(['gene_id']).sum()

df_genes_group_round = df_genes_group.apply(lambda x: x.round(0), 0).reset_index()

df_genes_group.to_csv('/nobackup/shilab/Data/novogene/edgeR_count_matrix.csv', index=True)
df_genes_group_round.to_csv('/nobackup/shilab/Data/novogene/edgeR_count_matrix_round.csv', index=False)

df2 = df.head().apply(lambda x: x.round(0), 0)

results = pd.read_table('/nobackup/shilab/Data/novogene/Counts/differentially_expressed_genes_edgeR_R697_IGFBP7_silencing_gencodeV24_significant_only_cpm_cut.txt').reset_index()

mart = pd.read_table('/nobackup/shilab/Data/novogene/Counts/gencode.v25.annotation.gtf').rename(columns={'Ensembl Gene ID':'index'})

results_mart = results.merge(mart, on='index')

results_mart.to_excel('/nobackup/shilab/Data/novogene/Counts/differentially_expressed_genes_edgeR_R697_IGFBP7_silencing_gencodeV24_significant_only_cpm_cut_named_phenotypes.xlsx', index=False)