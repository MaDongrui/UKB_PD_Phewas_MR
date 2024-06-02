#!/bin/bash


cd /home/mdr/PHESANT/WAS/data
python3

import pandas as pd

df_a = pd.read_csv('QC_data.csv')
df_b = pd.read_csv('raw_ukbdata.csv', escapechar='\\')
df_merged = pd.merge(df_a, df_b, on='userId', how='left')
df_merged.to_csv('ukbdata.csv', index=False)

exit()

parallel -j 20 'Rscript phenomeScan.r \
--phenofile="/home/mdr/PHESANT/WAS/data/ukbdata.csv" \
--traitofinterestfile="/home/mdr/PHESANT/WAS/data/StandardPRS_PD.csv" \
--variablelistfile="/home/mdr/PHESANT/variable-info/outcome-info.tsv" \
--datacodingfile="/home/mdr/PHESANT/variable-info/data-coding-ordinal-info.txt" \
--traitofinterest="exposure" \
--resDir="/home/mdr/PHESANT/WAS/results/" \
--userId="userId" \
--sensitivity \
--genetic=TRUE \
--partIdx={1} \
--numParts=20' ::: {1..20}


cd /home/mdr/PHESANT/resultsProcessing/

Rscript mainCombineResults.r \
--resDir="/home/mdr/PHESANT/WAS/results/" \
--variablelistfile="/home/mdr/PHESANT/variable-info/outcome-info.tsv" \
--numParts=20
