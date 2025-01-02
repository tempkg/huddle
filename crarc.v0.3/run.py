"""
  Hua Sun
  1/2/25 v0.3
"""


import argparse
import os
import re
import pandas as pd

# collect input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model', type=str, default='local', help='local/local-fc/lsf')
parser.add_argument('-ng', type=int, default=200, help='--min-gex-count number')
parser.add_argument('-na', type=int, default=8500, help='--min-atac-count number')
parser.add_argument('-t', '--table', type=str, default='', help='data table for run multiple samples')
parser.add_argument('-o', '--outdir', type=str, default='out_cellranger_arc', help='out dir')

args = parser.parse_args()




def Main():
    path = os.path.dirname(__file__)
    fconfig = f'{path}/config.ini'
    fbsub = f'{path}/hpc.bsub.sh'
    fscript = f'{path}/crarc.count.sh'

    CheckFile('hpc.bsub', fbsub)
    CheckFile('script', fscript)
    CheckFile('config', fconfig)

    if (args.outdir != '') & (args.outdir != '.'):
        if not os.path.isdir(args.outdir):
            os.mkdir(args.outdir)
    
    # cellranger-arc count
    CellrangerARC_Count_Multiple(fbsub, fscript, fconfig, args.model, args.ng, args.na, args.table, args.outdir)
    



"""
    Set Func.
"""

## Check file (absolute path)
def CheckFile(tag, f_path):
    if not os.path.isfile(f_path):
        print(f'[ERROR] The {tag} {f_path} does not exist !')
        exit()



## Cellranger-arc count for multiple samples
def CellrangerARC_Count_Multiple(fbsub, fscript, fconfig, model, min_gex_count, min_atac_count,ftable, outdir):
    info = pd.read_csv(ftable, sep='\t', names=['ref', 'name', 'library'])

    for index, row in info.iterrows():
        ref = row['ref']
        name = row['name']
        name = name.replace('.', '_')
        library = row['library']
        
        print(ref, name, library)

        cmd = f'sh {fbsub} 64 16 arc_{name} bash {fscript} -C {fconfig} -R {model} -G {ref} -L {library} -N {name} -g {min_gex_count} -a {min_atac_count} -O {outdir}'
        os.system(cmd)




if __name__ == '__main__':
    Main()



