# Hua Sun
# v2.2


import argparse
import pandas as pd
import os
import re


# collect input arguments
parser = argparse.ArgumentParser()

parser.add_argument('-s', '--sample', type=str, default='sample', help='sample name')
parser.add_argument('-b', '--barcode', type=str, default='datasetbarcode_extracted.fastq', help='path to extracted fasq files containing lineage barcodes')
parser.add_argument('-o', '--outdir', type=str, default='out_lineage', help='out dir')

args = parser.parse_args()




def Main():
    sample = args.sample
    outdir = args.outdir
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # input
    EXTRACTED_FASTQ = args.barcode

    # output files
    # v1
    REFORMATED_FASTA = f'{outdir}/dataset_barcode_reformat.fa'
    BARCODE_TXT = f'{outdir}/barcodes.txt'
    LINEAGE_BARCODE_TXT = f'{outdir}/lineage_barcodes.txt'

    # v2
    CELL_LINEAGE_BARCODE_ALL = f'{outdir}/cell_lineageBarcodes.txt'
    CELL_LINEAGE_BARCODE_FILTERED = f'{outdir}/cell_lineageBarcodes.filtered.txt'
    CELL_LINEAGE_BARCODE_FILTERED_CELL_WITH_ONE_LINEAGE_BARCODE = f'{outdir}/cell_lineageBarcodes.filtered.cell_withOneLineageBarcode.txt'


    #------- S1 raw
    # 1.1 make fa
    with open(EXTRACTED_FASTQ, 'r') as f:
        barcodeList = f.readlines()

    # make fasta file for cell + lineage
    for i in range(0, len(barcodeList), 4):
        barcodeList[i] = ">" + sample + ',' + FindBetween(barcodeList[i], "_", "_" ) + "," + FindBetween_R(barcodeList[i], "_", " ")


    with open(REFORMATED_FASTA, "w") as f:
        for i in range(0,len(barcodeList), 4):
            f.write(str(barcodeList[i]) + "\n" + str(barcodeList[i+1]))

    del barcodeList
    print('[INFO] Making fasta!')


    cmd = f'grep ^\\> {REFORMATED_FASTA} | sed \'s/>//\' | sort -u > {BARCODE_TXT}'
    os.system(cmd)


    cmd = f'grep -v \'^>\' {REFORMATED_FASTA} | sort | uniq -c | perl -pe \'s/^ +//\' | perl -pe \'s/ +/\\t/\' | sort -k1,1nr > {LINEAGE_BARCODE_TXT}'
    os.system(cmd)



    #------- S2 filtered
    with open(EXTRACTED_FASTQ, 'r') as f:
        fq_barcode = f.readlines()

    for i in range(0, len(fq_barcode), 4):
        fq_barcode[i] = FindBetween(fq_barcode[i], "_", "_" ) + "," + FindBetween_R(fq_barcode[i], "_", " ")
    
    with open(CELL_LINEAGE_BARCODE_ALL, "w") as f:
        for i in range(0, len(fq_barcode), 4):
            f.write(ExtractionCell(str(fq_barcode[i])) + "\t" + str(fq_barcode[i+1]))

    del fq_barcode
    print('[INFO] Created cell lineage barcodes!')


    lineage_barcode = pd.read_csv(CELL_LINEAGE_BARCODE_ALL, sep='\t', header=None)
    lineage_barcode.columns = ['cell', 'barcodes']

    lineage_barcode = lineage_barcode.drop_duplicates()
    lineage_barcode.to_csv(CELL_LINEAGE_BARCODE_ALL, sep='\t', header=None, index=False)

    lineage_barcode = lineage_barcode[lineage_barcode["barcodes"].str.contains('..CTG..ACT..GAC..TGA..CTG..ACT..GAC..')]

    lineage_barcode = lineage_barcode[~lineage_barcode["barcodes"].str.contains('N')]
    lineage_barcode.to_csv(CELL_LINEAGE_BARCODE_FILTERED, sep='\t', header=None, index=False)
    print('[INFO] Completed filtering!')


    cell_count = lineage_barcode.groupby(['cell']).count()
    cell_count = cell_count.reset_index()
    cell_count.columns = ['cell', 'Freq']
    
    cell_count = cell_count.loc[(cell_count['Freq'] == 1)]
    cell_list = list(cell_count.cell)

    lineage_barcode = lineage_barcode[lineage_barcode['cell'].isin(cell_list)]
    lineage_barcode.to_csv(CELL_LINEAGE_BARCODE_FILTERED_CELL_WITH_ONE_LINEAGE_BARCODE, sep='\t', header=None, index=False)

    print('[INFO] Completed !')




'''
    Set Func.
'''

## Find_between
def FindBetween( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""


## Find_between_R
def FindBetween_R( s, first, last ):
    try:
        start = s.rindex( first ) + len( first )
        end = s.rindex( last, start )
        return s[start:end]
    except ValueError:
        return ""


def ExtractionCell( s ): 
    return re.sub(',\\w+', '', s)





if __name__ == '__main__':
    Main()





