
import argparse
import pandas as pd
import os
import re


# collect input arguments
parser = argparse.ArgumentParser()

parser.add_argument('--cluster', action='store_true', default=False, help='use it in server')

parser.add_argument('-d', '--dir', type=str, default='.', help='processed dir')

parser.add_argument('-n', '--name', type=str, default='', help='folder name')
parser.add_argument('-o', '--outdir', type=str, default='out_lineage', help='out dir')

args = parser.parse_args()




def Main():
    cluster = args.cluster

    path = os.path.dirname(__file__)
    f_bsub = f'{path}/src/hpc.bsub.sh'
    src = f'{path}/src'
    script = f'{src}/reformatBarcodes.py'

    if (args.outdir != '') & (args.outdir != '.'):
        if not os.path.isdir(args.outdir):
            os.mkdir(args.outdir)
    
    # single sample
    if args.name != '':
        cmd = f'python3 {script} -s {args.name} -b {args.dir}/dataset_barcode_extracted.fastq -o {args.outdir}'
        # 'cluster' mode
        if cluster:
            cmd = f'sh {f_bsub} 526 1 rf_{args.name} python3 {script} -s {args.name} -b {args.dir}/dataset_barcode_extracted.fastq -o {args.outdir}'
            
        os.system(cmd)


    # multiple samples in one folder
    samples = [x for x in os.listdir(args.dir) if os.path.isdir(f'{args.dir}/{x}')]

    for x in samples:
        outdir_sample = f'{args.outdir}/{x}'
        cmd = f'python3 {script} -s {x} -b {args.dir}/{x}/dataset_barcode_extracted.fastq -o {outdir_sample}'
        # 'cluster' mode
        if cluster:
            cmd = f'sh {f_bsub} 526 1 rf_{x} python3 {script} -s {x} -b {args.dir}/{x}/dataset_barcode_extracted.fastq -o {outdir_sample}'
            
        os.system(cmd)




if __name__ == '__main__':
    Main()





