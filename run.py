"""
  Hua Sun
"""


import argparse
import yaml
import os


# collect input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--dir', type=str, default='', help='directory path of cellranger results for multiple samples')
parser.add_argument('--min_cells', type=int, default=3, help='cutoff min cells')
parser.add_argument('-r', '--ref', type=str, default='mm10', help='reference')
parser.add_argument('--yml', type=str, default='config.yml', help='yaml file')
parser.add_argument('-o', '--outdir', type=str, default='out_seurat_signac_rds', help='out dir')

args = parser.parse_args()




def Main():
    with open(args.yml, 'r') as f:
        data = yaml.safe_load(f)

    if not os.path.isdir(args.outdir):
        os.mkdir(args.outdir)

    RunSeuratSignacForMultiome(f_bsub=data['bsub'], r_scr=data['rscr'], run_scr=data['scr'], dir_arc=args.dir, ref=args.ref, min_cells=args.min_cells, outdir=args.outdir)


    



"""
    Set Func.
"""

def CheckFile(tag, f_path):
    if not os.path.isfile(f_path):
        print(f'[ERROR] The {tag} {f_path} does not exist !')
        exit()



def RunSeuratSignacForMultiome(f_bsub, r_scr, run_scr, dir_arc, ref, min_cells, outdir):
    sample_lists = [f.name for f in os.scandir(dir_arc) if f.is_dir()]

    for sample in sample_lists:
        h5 = f'{dir_arc}/{sample}/outs/filtered_feature_bc_matrix.h5'
        atac_fragment = f'{dir_arc}/{sample}/outs/atac_fragments.tsv.gz'
        barcode_metrics = f'{dir_arc}/{sample}/outs/per_barcode_metrics.csv'

        CheckFile('h5', h5)
        CheckFile('atac fragment', atac_fragment)
        CheckFile('barcode metrics', barcode_metrics)

        out_seurat = f'{outdir}/{sample}'
        if not os.path.isdir(out_seurat):
            os.mkdir(out_seurat)

        cmd = f'sh {f_bsub} 256 1 {sample} {r_scr} {run_scr} \
            --name {sample} --ref {ref} --h5 {h5} --min_cells {min_cells} \
            --atac_fragment {atac_fragment} \
            --barcode_metrics {barcode_metrics} \
            --outdir {out_seurat}'

        os.system(cmd)







if __name__ == '__main__':
    Main()



