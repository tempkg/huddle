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
parser.add_argument('-v', '--ver', type=str, default='v9', help='cellranger version v7.2/v9')
parser.add_argument('-m', '--model', type=str, default='local', help='local/local-fc/lsf')
parser.add_argument('-nc', type=int, default=10000, help='--force-cells number')
parser.add_argument('--chem', type=str, default='auto', help='--chemistry')
parser.add_argument('-t', '--table', type=str, default='input.txt', required=True, help='data table for run multiple samples')

parser.add_argument('-o', '--outdir', type=str, default='.', help='out dir')

args = parser.parse_args()




def Main():
    path = os.path.dirname(__file__)
    fconfig = f'{path}/config.{args.ver}.ini'
    src = f'{path}/src'
    fbsub = f'{path}/src/hpc.bsub.sh'

    CheckDir('src', src)
    CheckFile('config', f_config)

    outdir = args.outdir
    if (outdir == '.'):
        outdir = f'out_cellranger.{args.ver}'

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    

    # set script
    fscript = f'{src}/cellranger_{args.ver}.count.sh'
    CheckFile('script', fscript)

    CallSeqCount_MultipleSamples(fbsub, fconfig, fscript, args.model, args.nc, args.chem, args.table, outdir)    

    



"""
    Set Func.
"""

## Check file (absolute path)
def CheckFile(tag, f_path):
    if not os.path.isfile(f_path):
        print(f'[ERROR] The {tag} {f_path} does not exist !')
        exit()

## Check dir (absolute path)
def CheckDir(tag, dir_path):
    if not os.path.isdir(dir_path):
        print(f'[ERROR] The {tag} {dir_path} does not exist !')
        exit()



## CallSeqCount_MultipleSamples
def CallSeqCount_MultipleSamples(f_bsub, f_config, script, model, ncell, chemistry, f_table, outdir):

    info = pd.read_csv(f_table, sep='\t', names=['ref', 'name', 'fq_dir'])

    for index, row in info.iterrows():
        ref = row['ref']
        name = row['name']
        name = name.replace('.', '_')
        fq_dir = row['fq_dir']
        
        print(ref, name, fq_dir)

        cmd = f'sh {f_bsub} 64 16 scrna_{name} bash {script} -C {f_config} -R {model} -n {ncell} -c {chemistry} -G {ref} -F {fq_dir} -N {name} -O {outdir}'
        os.system(cmd)






if __name__ == '__main__':
    Main()



