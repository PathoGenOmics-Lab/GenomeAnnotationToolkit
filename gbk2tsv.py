import os
import glob
import pathlib
from Bio import SeqIO, SeqFeature
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(description = "Convert GenBank files to tab-delimited text files")
    parser.add_argument("-g", "--gbk", nargs = "+", type = str, required = True, dest = "gbks", default = "", help = "Input GenBank files")
    parser.add_argument("-o", "--outdir", type = str, required = False, dest = "outdir", default = ".", help = "Output directory")
    parser.add_argument("-f", "--features", type = str, required = False, dest = "features", default = "CDS,tRNA,rRNA", help = "Comma-separated features to store")
    parser.add_argument("-n", "--nucl_seq", action = "store_true", required = False, dest = "nucl_seq", help = "Turn on this option to print nucleotide sequences of features")
    parser.add_argument("-p", "--prot_seq", action = "store_true", required = False, dest = "prot_seq", help = "Turn on this option to print protein sequences of CDS")    
    return parser.parse_args()

def get_input_filenames(gbks):
    gbk_list = list(gbks)
    if len(gbk_list) == 1 and gbk_list[0].startswith("*"):  # *.gbk
        gbk_list = glob.glob(str(pathlib.Path(".", gbk_list[0])))  # Get names of all GenBank files under the current working directory
    return gbk_list

def main():
    args = parse_args()
    gbk_list = get_input_filenames(args.gbks)
    features = args.features.split(",")
    
if __name__ == "__main__":
    main()
