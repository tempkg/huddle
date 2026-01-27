"""
  Hua Sun
  v3.0
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
parser.add_argument('-y', '--yml', type=str, default='', help='yaml file')
parser.add_argument('-r', '--ref', type=str, default='mouse', help='reference')
parser.add_argument('-o', '--outdir', type=str, default='out_seurat5', help='out dir')

args = parser.parse_args()




def Main():
    path = os.path.dirname(__file__)

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    RunSeurat(
        f_bsub=data['bsub'], 
        hpc_mem=args.mem, 
        hpc_core=args.core, 
        r_scr=data['rscr'], 
        rcode=data['rcode'], 
        ftable=args.table, 
        dir_arc=args.dir, 
        ref=args.ref, 
        yml=args.yml,
        outdir=args.outdir
    )


"""
    Set Func.
"""

def CheckFile(tag, f_path):
    if not os.path.isfile(f_path):
        print(f'[ERROR] The {tag} {f_path} does not exist !')
        exit()



def RunSeurat(
        f_bsub, 
        hpc_mem, 
        hpc_core, 
        r_scr, 
        rcode, 
        ftable, 
        dir_arc, 
        ref, 
        yml,
        outdir='.'):
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

        # h5 file
        fh5 = f'{sample_dir}/filtered_feature_bc_matrix.h5'

        out_seurat = f'{outdir}/{sample}'

        cmd = f'sh {f_bsub} {hpc_mem} {hpc_core} seu_{sample} {r_scr} {rcode} --fh5 {fh5} \
            --name {sample} --ref {ref} --yaml {yml} --outdir {out_seurat}'

        os.system(cmd)







if __name__ == '__main__':
    Main()



