"""
  Hua Sun
  v2
"""


import argparse
import yaml
import pandas as pd
import os


# collect input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--mem', type=int, default=16, help='hpc memory')
parser.add_argument('--core', type=int, default=2, help='hpc core')
parser.add_argument('-d', '--dir', type=str, default='', help='directory path of cellranger results for multiple samples')
parser.add_argument('-t', '--table', type=str, default='', help='table')
parser.add_argument('--min_cells', type=int, default=10, help='cutoff min cells')
parser.add_argument('-r', '--ref', type=str, default='mouse', help='reference')
parser.add_argument('-o', '--outdir', type=str, default='out_signac', help='out dir')

args = parser.parse_args()




def Main():
    path = os.path.dirname(__file__)
    fyaml = f'{path}/config.yaml'

    with open(fyaml, 'r') as f:
        data = yaml.safe_load(f)

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    RunSignac(f_bsub=data['bsub'], hpc_mem=args.mem, hpc_core=args.core, r_scr=data['rscr'], run_scr=data['scr'], ftable=args.table, dir_arc=args.dir, ref=args.ref, min_cells=args.min_cells, outdir=args.outdir)


"""
    Set Func.
"""

def CheckFile(tag, f_path):
    if not os.path.isfile(f_path):
        print(f'[ERROR] The {tag} {f_path} does not exist !')
        exit()



def RunSignac(f_bsub, hpc_mem, hpc_core, r_scr, run_scr, ftable, dir_arc, ref, min_cells=10, outdir='.'):
    df = None
    sample_lists = None

    if dir_arc != '':
        sample_lists = [f.name for f in os.scandir(dir_arc) if f.is_dir()]
    else:
        df = pd.read_csv(ftable, sep='\t', header=None)
        sample_lists = df.iloc[:,0].tolist()


    for sample in sample_lists:
        if dir_arc != '':
            sample_dir = f'{dir_arc}/{sample}/outs'
        else:
            sample_dir = df.loc[df[0]==sample, 1].iloc[0]

        out_signac = f'{outdir}/{sample}'

        cmd = f'sh {f_bsub} {hpc_mem} {hpc_core} seu_{sample} {r_scr} {run_scr} --dir {sample_dir} \
            --name {sample} --ref {ref} --min_cells {min_cells} \
            --outdir {out_signac}'

        os.system(cmd)







if __name__ == '__main__':
    Main()



