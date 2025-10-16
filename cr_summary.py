"""
  Hua Sun
  4/29/25 v0.4
  
  QC for 10xGenomic cellranger

"""


import glob
import sys
import pandas as pd
import os
import re
import argparse



parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dir', type=str, default='out_cellranger', help='cellranger outs')
parser.add_argument('-o', '--outdir', type=str, default='out_cellranger_qc', help='out dir')

args = parser.parse_args()



def Main():
    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    # summary per sample
    samples = get_first_level_folders(args.dir)
    for x in samples:
        print(x)
        f_csv = f'{args.dir}/{x}/outs/metrics_summary.csv'
        QC_CellrangerSummaryFile(x, f_csv, f'{args.outdir}/{x}.metrics_summary.log')

    # merge
    MergeSummaryFile(args.outdir)



"""
  Set Func.
"""

def get_first_level_folders(path):
    folders = []
    for entry in os.scandir(path):
        if entry.is_dir():
            folders.append(entry.name)
    return folders



def Convert_to_numeric(value):
    value = str(value)

    if '%' in value:
        return round(float(value.replace('%', '')) / 100, 3)
    else:
        value = value.replace(',', '')
        
        if '.' not in value:
            return int(float(value))
        else:
            return float(value)



## Cellragner-arc count 'summary.csv'
def QC_CellrangerSummaryFile(name, f_data, outfile):
    df = pd.read_csv(f_data)
    df = df.T
    df.columns = [name]
    print(df[name].tolist())
    df[name] = [Convert_to_numeric(x) for x in df[name].tolist()]
    #print(df)
    report = pd.DataFrame(
        [['Estimated Number of Cells','500-10000','','WARNING'],
        ['Mean Reads per Cell','>=20000','','WARNING'],
        ['Median Genes per Cell','>1000','','WARNING'],
        ['Number of Reads','-','','-'],
        ['Valid Barcodes','>0.75','','WARNING'],
        ['Valid UMIs','>0.75','','WARNING'],
        ['Sequencing Saturation','-','','-'],
        ['Q30 Bases in Barcode','>0.85','','WARNING'],
        ['Q30 Bases in RNA Read','>0.85','','WARNING'],
        ['Q30 Bases in UMI','>0.85','','WARNING'],
        ['Reads Mapped to Genome','>0.85','','WARNING'],
        ['Reads Mapped Confidently to Genome','>0.5','','FAIL'],
        ['Reads Mapped Confidently to Intergenic Regions','<0.3','','WARNING'],
        ['Reads Mapped Confidently to Intronic Regions','<0.3','','FAIL'],
        ['Reads Mapped Confidently to Exonic Regions','>0.3','','FAIL'],
        ['Reads Mapped Confidently to Transcriptome','>0.5','','WARNING'],
        ['Reads Mapped Antisense to Gene','<0.1','','WARNING'],
        ['Fraction Reads in Cells','>0.7','','WARNING'],
        ['Total Genes Detected','-','','-'],
        ['Median UMI Counts per Cell','100','','WARNING']],

        columns=['Catalog', '10xStandard', name, f'{name}.qc'])

    for catalog in report.Catalog:
      val = df.loc[catalog, name]
      report.loc[report.Catalog==catalog, name] = df.loc[catalog, name]

      if catalog == 'Estimated Number of Cells':
          if (val > 500) & (val < 10000):
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 500:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Mean Reads per Cell':
          if val >= 20000:
              print(val)
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
      if catalog == 'Median Genes per Cell':
          if val > 1000:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 200:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Valid Barcodes':
          if val > 0.75:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 0.5:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Valid UMIs':
          if val > 0.75:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 0.5:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Q30 Bases in Barcode':
          if val > 0.85:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 0.3:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Q30 Bases in RNA Read':
          if val > 0.85:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 0.3:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Q30 Bases in UMI':
          if val > 0.85:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 0.3:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Reads Mapped to Genome':
          if val > 0.85:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 0.5:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Reads Mapped Confidently to Genome':
          if val > 0.5:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS' 
          if val < 0.3:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'    
      if catalog == 'Reads Mapped Confidently to Intergenic Regions':
          if val < 0.3:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
      if catalog == 'Reads Mapped Confidently to Intronic Regions':
          if val < 0.3:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
      if catalog == 'Reads Mapped Confidently to Exonic Regions':
          if val > 0.3:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
      if catalog == 'Reads Mapped Confidently to Transcriptome':
          if val > 0.5:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 0.3:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Reads Mapped Antisense to Gene':
          if val < 0.1:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'  
          if val > 0.4:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'    
      if catalog == 'Fraction Reads in Cells':
          if val > 0.7:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'
          if val < 0.4:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'FAIL'
      if catalog == 'Median UMI Counts per Cell':
          if val > 100:
              report.loc[report.Catalog==catalog,f'{name}.qc'] = 'PASS'

    #print(report)
    n_pass = report.loc[report.iloc[:,3]=='PASS'].shape[0]/report.loc[report.iloc[:,3].isin(['PASS','WARNING','FAIL'])].shape[0] * 100
    n_pass = round(n_pass, 2)
    print(f'Percentage of PASS: {n_pass}%')
    report = report._append(dict(zip(report.columns,['Percentage of PASS', '.', '.', n_pass])), ignore_index=True)
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









