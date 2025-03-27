"""
    Hua Sun
    v0.1

"""

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument('-r', '--r_script', default='~/software/miniconda3/envs/r442/bin/Rscript', help='Rscript location')
parser.add_argument('-o', '--outdir', default='out_infercnv_sec', help='out dir')
args = parser.parse_args()


def main():
    path = os.path.dirname(__file__)
    r_script = args.r_script

    print('[INFO] Call CNV ...')
    fscript = f'{path}/re.infercnv.R'
    fconfig = f'{path}/config.yml'
    cmd = f'{r_script} {fscript} {fconfig} {args.outdir}'
    os.system(cmd)



if __name__ == '__main__':
    main()


