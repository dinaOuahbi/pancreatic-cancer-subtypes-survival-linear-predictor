# DO
#31-01-2023

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#LOADING
OUT = ''
groups = pd.read_csv('groups.csv') # PATIENT WITH HIGH OR LOW
workbook = pd.read_excel("paad_subtype.xlsx") # SUBTYPE DATA FROM GDC PORTAL
lin_pred = pd.read_csv('condition.csv', index_col=0)


# PROCESS SUBTYPE DATA
df = workbook[['Tumor Sample ID','mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX']]
df = df.rename(columns={'Tumor Sample ID':'patient',
                  'mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX':'cancer_subtype'})
df['patient']=['-'.join(i.split('-')[:3]).strip() for i in df['patient']]
print(f'shape OF SUBTYPE : {groups.shape}')

# first merge with high low groups ===> fisher
merge = pd.merge(groups, df, on='patient')
merge.shape
print(pd.crosstab(merge['pred_dicho'],merge['cancer_subtype'],dropna=False))


# predicteur lineaire
lin_pred = lin_pred[['patient','lbl_pred_all']]
lin_pred['patient']=[i[:12].strip() for i in lin_pred['patient']] # if MIL model
lin_pred.drop_duplicates('patient', inplace=True)
print(f'shape of linear pred : {lin_pred.shape}')


# seconde merge with pred lineaire ==> boxplot
merge2 = pd.merge(lin_pred, merge, on='patient')
print(f'shape of merge2 : {merge2.shape}')

merge2.to_csv(f'{OUT}paad_subtype.csv')
