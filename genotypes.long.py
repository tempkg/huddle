import pandas as pd
import argparse

# collect input arguments
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--name', type=str, default='sample', help='name')
parser.add_argument('-i', '--input', type=str, default='genotype.tsv', help='genotype file')
parser.add_argument('-o', '--output', type=str, default='file.out', help='out file')

args = parser.parse_args()



def convert_hla_format_pandas(name="sample1", input_file='in.tsv', output_file='out.txt'):
    """
    Convert HLA data from wide format to long format using pandas
    """
    
    # Read the data
    df = pd.read_csv(input_file, sep='\s+')
    
    # Reshape from wide to long format
    df_long = df.melt(id_vars=['subject'], 
                      var_name='gene', 
                      value_name='hla')
    
    # Remove 'subject' column and rename
    df_long = df_long.drop('subject', axis=1)
    df_long.insert(0, 'sample', name)
    
    # Save to file
    df_long.to_csv(output_file, sep='\t', index=False)


# Usage
convert_hla_format_pandas(args.name, args.input, args.output)

