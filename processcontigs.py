import os
import argparse
from collections import defaultdict

def parse_contigs(input_file, output_folder):
    # Dictionary to hold the contigs for each unique assembly
    contigs_dict = defaultdict(list)
    
    # Read the input file and group lines by assembly
    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                # Extract the unique assembly identifier (donor#hap)
                assembly = '#'.join(line.split('#')[:2])
                contigs_dict[assembly].append(line)
    
    # Write the grouped lines to separate files
    for assembly, contigs in contigs_dict.items():
        output_file = os.path.join(output_folder, f'chr1.{assembly}.contig')
        with open(output_file, 'w') as out_file:
            for contig in contigs:
                out_file.write(contig + '\n')

def main():
    parser = argparse.ArgumentParser(description="Split contigs file into separate assembly files")
    parser.add_argument('-I', '--input', required=True, help="Input contigs file")
    parser.add_argument('-o', '--output', required=True, help="Output folder for separate assembly files")
    
    args = parser.parse_args()
    
    # Create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    
    # Parse the contigs file and create separate files
    parse_contigs(args.input, args.output)

if __name__ == "__main__":
    main()
