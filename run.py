"""
    Hua Sun
    2024-11-18 v0.32

"""

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, default='', help='input file')
parser.add_argument('-r', '--r_script', default='~/software/miniconda3/envs/r442/bin/Rscript', help='Rscript location')
args = parser.parse_args()


def main():
    path = os.path.dirname(__file__)
    r_script = args.r_script

    # 1
    print('[INFO] Extraction data ...')
    fscript = f'{path}/src/seu2mtx.R'
    cmd = f'{r_script} {fscript} {args.input}'
    os.system(cmd)

    # 2
    print('[INFO] Call CNV ...')
    fscript = f'{path}/src/infercnv.R'
    cmd = f'{r_script} {fscript} {args.input}'
    os.system(cmd)



if __name__ == '__main__':
    main()


