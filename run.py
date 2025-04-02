"""
  Hua Sun
"""


import argparse
import yaml
import os


# collect input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--dir', type=str, default='', help='directory path of cellranger results for multiple samples')
parser.add_argument('--min_cells', type=int, default=3, help='cutoff min cells')
parser.add_argument('-r', '--ref', type=str, default='mm10', help='reference')
parser.add_argument('-o', '--outdir', type=str, default='out_seurat5', help='out dir')

args = parser.parse_args()




def Main():
    path = os.path.dirname(__file__)
    fyaml = f'{path}/config.yaml'

    with open(fyaml, 'r') as f:
        data = yaml.safe_load(f)

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    RunSeurat(f_bsub=data['bsub'], r_scr=data['rscr'], run_scr=data['scr'], dir_arc=args.dir, ref=args.ref, min_cells=args.min_cells, outdir=args.outdir)


"""
    Set Func.
"""

def CheckFile(tag, f_path):
    if not os.path.isfile(f_path):
        print(f'[ERROR] The {tag} {f_path} does not exist !')
        exit()



def RunSeurat(f_bsub, r_scr, run_scr, dir_arc, ref, min_cells, outdir):
    sample_lists = [f.name for f in os.scandir(dir_arc) if f.is_dir()]

    for sample in sample_lists:
        h5 = f'{dir_arc}/{sample}/outs/filtered_feature_bc_matrix.h5'

        CheckFile('h5', h5)

        out_seurat = f'{outdir}/{sample}'

        cmd = f'sh {f_bsub} 10 2 {sample} {r_scr} {run_scr} \
            --name {sample} --ref {ref} --h5 {h5} --min_cells {min_cells} \
            --outdir {out_seurat}'

        os.system(cmd)







if __name__ == '__main__':
    Main()



