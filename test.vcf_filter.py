

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', type=str, default='sample1.vcf', help='cellsnp vcf')
parser.add_argument('-o', '--out', type=str, default='formated.tsv', help='output file')
parser.add_argument('-n', '--name', type=str, default='sample', help='sample name')
parser.add_argument('--AD', type=int, default=3, help='min AD')
parser.add_argument('--DP', type=int, default=10, help='min DP')

args = parser.parse_args()



vcf = args.input
out = args.out

name = args.name
minAD = args.AD
minDP = args.DP


with open(vcf, 'r') as vcf, open(out, 'w') as tsv:
    for line in vcf:
        if line.startswith('#'):
            continue

        parts = line.strip().split('\t')
        chrom, pos, _, ref, alt, _, _, info = parts
        
        ad = None
        dp = None
        for item in info.split(';'):
            if item.startswith('AD='):
                ad = int(item.split('=')[1])
            elif item.startswith('DP='):
                dp = int(item.split('=')[1])
        
        if ad >= minAD and dp >= minDP:
            vaf = round((ad / dp) * 100, 2)
            
            if vaf == int(vaf):
                vaf = str(int(vaf))
            else:
                vaf = str(vaf)
            
            tsv.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{vaf}\t{name}\n")

print("Done! Check formated.tsv")




