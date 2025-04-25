"""
  Hua Sun
  9/25/24 v0.3
  
  QC for 10xGenomic single cell data

"""


import glob
import sys
import pandas as pd
import os
import re
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dir', type=str, default='.', help='in dir')
parser.add_argument('-o', '--outdir', type=str, default='out_cratac_qc', help='out dir')

args = parser.parse_args()



def Main():
  if not os.path.isdir(args.outdir):
      os.mkdir(args.outdir)

  f_csv = f'{args.dir}/outs/summary.csv'
  outfile = f'{args.outdir}/cratac_qc.summary.log'

  QC_CellrangerATACSummaryFile(f_csv, outfile)


  # summary per sample
  samples = get_first_level_folders(args.dir)
  for x in samples:
      print(x)
      f_csv = f'{args.dir}/{x}/outs/summary.csv'
      QC_CellrangerATACSummaryFile(x, f_csv, f'{args.outdir}/{x}.metrics_summary.log')
  

  # merge
  MergeSummaryFile(args.outdir)


"""
  Set Func.
"""

## Cellragner-arc count 'summary.csv'
def QC_CellrangerATACSummaryFile(f_data, outfile):
  df = pd.read_csv(f_data, header=None)
  df = df.T
  df.columns = df.iloc[0]
  df.drop(df.index[0], inplace=True)
  df.set_index('Sample ID', inplace=True)

  name = df.columns[0]
  report = pd.DataFrame(
      [['Estimated number of cells','500-10000','','FAIL'],
        ['Confidently mapped read pairs','>0.8','','FAIL'],
        ['Estimated bulk library complexity','-','','FAIL'],
        ['Fraction of all fragments in cells','-','','FAIL'],
        ['Fraction of all fragments that pass all filters and overlap called peaks','-','','FAIL'],
        ['Fraction of genome in peaks','<0.75','','FAIL'],
        ['Fraction of high-quality fragments in cells','>0.4','','FAIL'],
        ['Fraction of high-quality fragments overlapping TSS','-','','FAIL'],
        ['Fraction of high-quality fragments overlapping peaks','>0.25','','FAIL'],
        ['Fraction of transposition events in peaks in cells','>0.25','','FAIL'],
        ['Fragments in necleosome-free regions','-','','FAIL'],

        ['Mean raw read pairs per cell','>5000','','FAIL'],
        ['Median high-quality fragments per cell','>100','','FAIL'],
        ['Non-nuclear read pairs','<0.1','','FAIL'],
        
        ['Number of peaks','-','','-'],
        ['Percent duplicates','-','','-'],
        ['Q30 bases in barcode','-','','-'],
        ['Q30 bases in read 1','-','','-'],
        ['Q30 bases in read 2','-','','-'],
        ['Q30 bases in sample index i1','-','','-'],
        ['Sequenced read pairs','-','','-'],
        ['TSS enrichment score','>5','','FAIL'],
        ['Unmapped read pairs','-','','-'],
        ['Valid barcodes','>0.85','','FAIL']],

      columns=['Catalog', '10xStandard', name, f'{name}.qc'])

  for catalog in report.Catalog:
    val = df.loc[catalog, name]
    report.loc[report.Catalog==catalog,name] = df.loc[catalog, name]

    val = float(val)

    if catalog == 'Estimated number of cells':
      if (val > 500) & (val < 10000):
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'Confidently mapped read pairs':
      if val > 0.8:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'Fraction of genome in peaks':
      if val < 0.75:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'Fraction of high-quality fragments in cells':
      if val > 0.4:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'Fraction of high-quality fragments overlapping peaks':
      if val > 0.25:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'Fraction of transposition events in peaks in cells':
      if val > 0.25:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'Mean raw read pairs per cell':
      if val > 5000:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'Median high-quality fragments per cell':
      if val > 100:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'Non-nuclear read pairs':
      if val < 0.1:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'TSS enrichment score':
      if val > 5:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
    if catalog == 'Valid barcodes':
      if val > 0.85:
        report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'

  #print(report)
  #n_pass = report.loc[report.iloc[:,3]=='PASS'].shape[0]/report.shape[0] * 100
  n_pass = report.loc[report.iloc[:,3]=='PASS'].shape[0]/report.loc[report.iloc[:,3].isin(['PASS','FAIL'])].shape[0] * 100
  n_pass = round(n_pass, 2)
  print(f'Percentage of PASS: {n_pass}%')
  report = report.append(dict(zip(report.columns,['Percentage of PASS', '.', '.', n_pass])), ignore_index=True)
  report.to_csv(outfile, sep='\t', index=False)






# merge summary files
def MergeSummaryFile(dir_qc):
    files = glob.glob(f'{dir_qc}/*.metrics_summary.log')
    i = 0
    df_qc = []
    for f in files:
        d_temp = []
        if i == 0:
            d_temp = pd.read_csv(f, sep='\t')
            df_qc = d_temp
            i = 1
        else:
            d_temp = pd.read_csv(f, sep='\t', usecols=[0,2,3])
            df_qc = df_qc.merge(d_temp, on='Catalog')

        print(f'[INFO] file:{f} - {d_temp.shape}')

    df_qc.to_csv(f'{dir_qc}/merged_all.metrics_summary.xls', sep='\t', index=False)






if __name__ == '__main__':
  Main()






