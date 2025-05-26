import pandas as pd
import warnings
import argparse
import os

# Suppress warnings
warnings.filterwarnings("ignore")

def process_paf(input_file, output_file):
    # Load the PAF file into a DataFrame
    filteredpaf = pd.read_csv(input_file, sep='\t', header=None)
    originalpaf = filteredpaf.copy()

    # Split the column and extract the last part after the '#'
    filteredpaf[14] = filteredpaf[5].str.split('#').str[-1]

    # Filter to rows where target is duplicated
    filteredpaf = filteredpaf[filteredpaf.duplicated(subset=filteredpaf.columns[0], keep=False)]

    # Remove the duplicate tuple rows of "(target,chrlabel)"
    filteredpaf = filteredpaf.drop_duplicates(subset=[filteredpaf.columns[0], filteredpaf.columns[14]])

    # Keep only the tuples where target is duplicated with different chr labels
    filteredpaf = filteredpaf[filteredpaf.duplicated(subset=filteredpaf.columns[0], keep=False)]

    # Select one random alignment for each duplicated alignment
    randomalignments = filteredpaf.groupby(filteredpaf.columns[0]).apply(lambda x: x.sample(1)).reset_index(drop=True)

    # Drop temporary column used for chrlabel comparison
    randomalignments = randomalignments.drop(columns=[14])
    filteredpaf = filteredpaf.drop(columns=[14])

    # Subtract the duplicates with different chr labels from the original paf
    originalpaf_minus_filtered = pd.concat([originalpaf, filteredpaf]).drop_duplicates(keep=False)

    # Concatenate the random alignments back into the originalpaf_minus_filtered
    originalpaf_added_random_alignments = pd.concat([originalpaf_minus_filtered, randomalignments], ignore_index=True)

    # Save the result to the output file
    originalpaf_added_random_alignments.to_csv(output_file, sep='\t', header=None, index=False)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Process PAF file to filter and randomize alignments.")
    parser.add_argument("-i", "--input", required=True, help="Input PAF file path")
    parser.add_argument("-o", "--output", required=True, help="Output PAF file path")

    # Parse arguments
    args = parser.parse_args()

    # Check if input file exists
    if not os.path.isfile(args.input):
        print(f"Error: The file {args.input} does not exist.")
        return

    # Process the input PAF file and save the output
    process_paf(args.input, args.output)
    print(f"Processing complete. Output saved to {args.output}")

if __name__ == "__main__":
    main()
