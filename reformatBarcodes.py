# Hua Sun
# v2.3


import argparse
import pandas as pd
import os
import re
import subprocess


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

    SUMMARY_LOG = f'{outdir}/summary.log'


    #------- S1 raw
    # 1.1 make fa
    with open(EXTRACTED_FASTQ, 'r') as f:
        barcodeList = f.readlines()

    # make fasta file for cell + lineage
    for i in range(0, len(barcodeList), 4):
        barcodeList[i] = ">" + sample + ',' + FindBetween(barcodeList[i], "_", "_" ) + "," + FindBetween_R(barcodeList[i], "_", " ")

    # output formated fasta - sample, barcode, umi with lineage barcode
    with open(REFORMATED_FASTA, "w") as f:
        for i in range(0,len(barcodeList), 4):
            f.write(str(barcodeList[i]) + "\n" + str(barcodeList[i+1]))

    del barcodeList
    print('[INFO] Making fasta!')


    # 1.2 output - sample, barcode, umi
    cmd = f'grep ^\\> {REFORMATED_FASTA} | sed \'s/>//\' | sort -u > {BARCODE_TXT}'
    os.system(cmd)
    print('[INFO] Extraction barcodes!')


    # 1.3 output - lineage & counts(reads)
    # show all without filter (but read counts > 9 may be trustable)
    cmd = f'grep -v \'^>\' {REFORMATED_FASTA} | sort | uniq -c | perl -pe \'s/^ +//\' | perl -pe \'s/ +/\\t/\' | sort -k1,1nr > {LINEAGE_BARCODE_TXT}'
    os.system(cmd)
    print('[INFO] Extraction lineage barcodes!')



    #------- S2 filtered
    # 2.1 make cell + lineage barcodes
    with open(EXTRACTED_FASTQ, 'r') as f:
        fq_barcode = f.readlines()

    for i in range(0, len(fq_barcode), 4):
        fq_barcode[i] = FindBetween(fq_barcode[i], "_", "_" ) + "," + FindBetween_R(fq_barcode[i], "_", " ")
    
    with open(CELL_LINEAGE_BARCODE_ALL, "w") as f:
        for i in range(0, len(fq_barcode), 4):
            f.write(ExtractionCell(str(fq_barcode[i])) + "\t" + str(fq_barcode[i+1]))

    del fq_barcode
    print('[INFO] Created cell lineage barcodes!')


    # 2.2 filter barcodes - only extract 100% match with consensus barcode
    lineage_barcode = pd.read_csv(CELL_LINEAGE_BARCODE_ALL, sep='\t', header=None)
    lineage_barcode.columns = ['cell', 'barcodes']
    # remove duplicated rows and overwite 'CELL_LINEAGE_BARCODE_ALL'
    # IMPORTANT - unique row 
    lineage_barcode = lineage_barcode.drop_duplicates()
    lineage_barcode.to_csv(CELL_LINEAGE_BARCODE_ALL, sep='\t', header=None, index=False)

    # only including conserved lineage barcodes - 37bp
    lineage_barcode = lineage_barcode[lineage_barcode["barcodes"].str.contains('..CTG..ACT..GAC..TGA..CTG..ACT..GAC..')]
    # filter barcodes with 'N'
    lineage_barcode = lineage_barcode[~lineage_barcode["barcodes"].str.contains('N')]
    lineage_barcode.to_csv(CELL_LINEAGE_BARCODE_FILTERED, sep='\t', header=None, index=False)
    print('[INFO] Completed filtering!')


    # 2.3 find cell with one lineage barcode
    # count frequency
    cell_count = lineage_barcode.groupby(['cell']).count()
    cell_count = cell_count.reset_index()
    cell_count.columns = ['cell', 'Freq']
    
    # cell with one lineage barcode
    cell_count = cell_count.loc[(cell_count['Freq'] == 1)]
    cell_list = list(cell_count.cell)

    # matched target cells
    lineage_barcode = lineage_barcode[lineage_barcode['cell'].isin(cell_list)]
    lineage_barcode.to_csv(CELL_LINEAGE_BARCODE_FILTERED_CELL_WITH_ONE_LINEAGE_BARCODE, sep='\t', header=None, index=False)

    print('[INFO] Completed !')


    # 2.4 Summary
    n_barcodes = CountRowsFromFile(LINEAGE_BARCODE_TXT)
    n_cells_all = CountUniqueCellsFromFile(CELL_LINEAGE_BARCODE_ALL)
    n_cells_filtered = CountUniqueCellsFromFile(CELL_LINEAGE_BARCODE_FILTERED)
    n_cells_with_one_barcode = CountUniqueCellsFromFile(CELL_LINEAGE_BARCODE_FILTERED_CELL_WITH_ONE_LINEAGE_BARCODE)

    data = {
        'Count': [n_barcodes, n_cells_all, n_cells_filtered, n_cells_with_one_barcode],
        'Catalog': ['All lineage-barcode', 'Cells with lineage-barcode', 'Cells with filtered lineage-barcode', 'Cells have only one filtered lineage-barcode']
    }

    df = pd.DataFrame(data)
    df.to_csv(SUMMARY_LOG, sep='\t', index=False)
    




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


## Extract cells
def ExtractionCell( s ): 
    return re.sub(',\\w+', '', s)


## CountRowsFromFile
def CountRowsFromFile(file):
    try:
        output = subprocess.check_output(
            ['wc', '-l', file],
            text=True
        )
        #print(f"Line count: {int(output.split()[0])}")
        line_count = int(output.split()[0])
        return line_count
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e}")


## CountUniqueCellsFromFile
def CountUniqueCellsFromFile(file):
    try:
        # Run the piped command safely
        result = subprocess.run(
            [f'cut -f1 {file} | sort -u | wc -l'],
            shell=True,
            capture_output=True,
            text=True,
            check=True
        )
        unique_count = int(result.stdout.strip())  # Extract the number
        return unique_count
    except subprocess.CalledProcessError as e:
        print(f"Command failed (exit {e.returncode}): {e.stderr}")
        return None







if __name__ == '__main__':
    Main()





