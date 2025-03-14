"""
    Hua Sun
    2024-11-18 v0.32

"""

import argparse
import pandas as pd
import re
import os
import yaml


# collect input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', required=True, default='', help='input file')
parser.add_argument('-l', '--loc', default='hpc', help='reference location')
args = parser.parse_args()


def main():
    loc = args.loc
    path = os.path.dirname(__file__)
    fconfig = f'{path}/src/config.yaml'
    config = ReadYaml(fconfig)
    r_script = 'Rscript'
    if 'RSCRIPT' in config[loc]:
        r_script = config[loc]['RSCRIPT']

    dy = ReadYaml(args.input)
    ref = dy['ref']
    outdir = dy['outdir']

    fgene = ''
    if (ref == 'mm10'):
        fgene = config[loc]['MM10_GENES']
    if (ref == 'hg38'):
        fgene = config[loc]['HG38_GENES']

    print(fgene)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)



    # 1
    print('[INFO] Extraction data ...')
    fscript = f'{path}/src/seu2mtx.R'
    cmd = f'{r_script} {fscript} {args.input}'
    os.system(cmd)

    # 2
    print('[INFO] Reformat data ...')
    ConvertData(fgene, outdir)

    # 3
    print('[INFO] Call CNV ...')
    fscript = f'{path}/src/infercnv.R'
    cmd = f'{r_script} {fscript} {args.input}'
    os.system(cmd)



'''
    subset
'''
def ReadYaml(f):
    with open(f, 'r') as file:
        data = yaml.safe_load(file)

    return data



def ConvertData(fgene, outdir):
    # read gene info
    gene_mtx = pd.read_csv(fgene, sep='\t', names=['gene', 'chr', 'start', 'end'])
    gene_mtx = gene_mtx.drop_duplicates(subset='gene', keep="first")

    # read cell info
    cell_anno = pd.read_csv(f'{outdir}/cell.anno', sep='\t', names=['cell', 'anno'])
    cells = cell_anno['cell'].tolist()
    
    # read count matrix (fast)
    print('[INFO] Reading count ...')
    chunk = pd.read_csv(f'{outdir}/count.mtx', sep='\t', header=None, chunksize=1000)
    tumor_mtx = pd.concat(chunk, ignore_index=True)
    with open(f'{outdir}/genes.txt') as f:
        rownames = [ line.strip() for line in f ]
    with open(f'{outdir}/barcodes.txt') as f:
        colnames = [ line.strip() for line in f ]
    tumor_mtx.columns = colnames
    tumor_mtx.index = rownames

    # extract
    print('[INFO] Extract ...')
    tumor_mtx = tumor_mtx.filter(items=gene_mtx.gene, axis=0)
    tumor_mtx = tumor_mtx.filter(items=cells, axis=1)
    print(f'[INFO] Tumor data after filter - {tumor_mtx.shape}')
    
    # re-extract gene & cells by matrix
    gene_mtx = gene_mtx.loc[gene_mtx['gene'].isin(tumor_mtx.index)]
    cell_anno = cell_anno.loc[cell_anno['cell'].isin(tumor_mtx.columns)]

    # output
    print('[INFO] Output ...')
    os.mkdir(f'{outdir}/mtx')
    tumor_mtx.to_csv(f'{outdir}/mtx/expCount.txt', sep='\t')
    cell_anno.to_csv(f'{outdir}/mtx/cellInfo.txt', sep='\t', index=False, header=None)
    gene_mtx.to_csv(f'{outdir}/mtx/geneLoci.txt', sep='\t', index=False, header=None)

    # remove processing data
    os.system(f'rm -f {outdir}/count.mtx')
    os.system(f'rm -f {outdir}/genes.txt')
    os.system(f'rm -f {outdir}/barcodes.txt')
    os.system(f'rm -f {outdir}/cell.anno')





if __name__ == '__main__':
    main()


