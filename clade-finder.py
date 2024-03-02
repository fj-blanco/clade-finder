"""
@author: fj-blanco
"""

import os
import argparse
from Common.CommonMethods import getJSON, getEncodedPositivesNegatives

def parse_args():
    parser = argparse.ArgumentParser(description='Process SNP and clade data to generate encoded positives and negatives, and output JSON data.')
    parser.add_argument('inputFileType',
                        type=str,
                        choices=['vcf', 'csv', 'txt'],
                        help='Type of the input file containing SNPs: vcf, csv, or txt. For txt, enter SNPs in a comma-separated format.'
                        )
    parser.add_argument('inputFilePath',
                        type=str,
                        help='Path to the input file containing SNP and/or Clade information.'
                        )
    parser.add_argument('cladeSNPFilePath',
                        type=str,
                        help='Path to the TB Clade SNP file.'
                        )
    parser.add_argument('SNPcladeFilePath',
                        type=str,
                        help='Path to the TB SNP clade file.'
                        )
    parser.add_argument('--params',
                        type=str,
                        default="score",
                        help=('Additional comma separated parameters for JSON generation, such as "all"'
                              'to include all possible clade predictions, not just the top one.')
                        )
    
    args = parser.parse_args()

    # Checking if the provided files exist
    files_to_check = [args.inputFilePath, args.cladeSNPFilePath, args.SNPcladeFilePath]
    for file_path in files_to_check:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file {file_path} does not exist.")
    
    return args

def main(input_type, input_file, clade_snp_file, snp_clade_file, params):
    if input_type == "txt":
        # read the input file and get the SNPs
        with open(input_file, "r") as f:
            snps = f.read().strip().split(",")

        (positives, negatives) = getEncodedPositivesNegatives(snps)

        # get the JSON CladeFinder's output
        output = getJSON(params, positives, negatives, clade_snp_file, snp_clade_file, None)
        print("---------------------")
        print("CladeFinder's prediction:", output)
        print("---------------------")
        return output
    else:
        raise ValueError("Input file type not supported yet.")

if __name__ == "__main__":
    args = parse_args()

    print("---------------------")
    print("Running CladeFinder with the following arguments:")
    print(f"Input type: {args.inputFileType}")
    print(f"Input file: {args.inputFilePath}")
    print(f"Clade SNP file: {args.cladeSNPFilePath}")
    print(f"SNP clade file: {args.SNPcladeFilePath}")
    print("---------------------")

    main(args.inputFileType, args.inputFilePath, args.cladeSNPFilePath, args.SNPcladeFilePath, args.params)