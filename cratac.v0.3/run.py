"""

  Hua Sun
  2024-12-20 v0.2.1

"""


import argparse
import os
import re
import pandas as pd

# collect input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-m', '--model', type=str, default='local', help='local/local-fc/lsf')
parser.add_argument('-nc', type=int, default=10000, help='--force-cells number')
parser.add_argument('--chem', type=str, default='auto', help='--chemistry')
parser.add_argument('-t', '--table', type=str, default='input.txt', required=True, help='data table for run multiple samples')

parser.add_argument('-o', '--outdir', type=str, default='out_cellranger_atac', help='out dir')

args = parser.parse_args()




def Main():
    path = os.path.dirname(__file__)
    fconfig = f'{path}/config.ini'
    fbsub = f'{path}/hpc.bsub.sh'
    fscript = f'{path}/cratac.count.sh'

    CheckFile('hpc.bsub', fbsub)
    CheckFile('script', fscript)
    CheckFile('config', fconfig)

    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    

    CallSeqCount_MultipleSamples(f_bsub, f_config, script, args.model, args.nc, args.chem, args.table, outdir)     

    



"""
    Set Func.
"""

def CheckFile(tag, f_path):
    if not os.path.isfile(f_path):
        print(f'[ERROR] The {tag} {f_path} does not exist !')
        exit()



def CallSeqCount_MultipleSamples(f_bsub, f_config, script, model, ncell, chemistry, f_table, outdir):

    info = pd.read_csv(f_table, sep='\t', names=['ref', 'name', 'fq_dir'])

    for index, row in info.iterrows():
        ref = row['ref']
        name = row['name']
        name = name.replace('.', '_')
        fq_dir = row['fq_dir']
        
        print(ref, name, fq_dir)

        cmd = f'sh {f_bsub} 64 16 snatac_{name} bash {script} -C {f_config} -R {model} -n {ncell} -c {chemistry} -G {ref} -F {fq_dir} -N {name} -O {outdir}'
        os.system(cmd)






if __name__ == '__main__':
    Main()



