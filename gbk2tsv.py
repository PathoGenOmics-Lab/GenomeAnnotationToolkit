import os
import glob
import pathlib
from Bio import SeqIO, SeqFeature
from argparse import ArgumentParser

def parse_args():
    parser = ArgumentParser(description = "Transforms GenBank files into tab-separated value files")
    parser.add_argument("-g", "--gbk", nargs = "+", type = str, required = True, dest = "gbks", default = "", help = "Specify GenBank files as input")
    parser.add_argument("-o", "--outdir", type = str, required = False, dest = "outdir", default = ".", help = "Define the directory for output files")
    parser.add_argument("-f", "--features", type = str, required = False, dest = "features", default = "CDS,tRNA,rRNA", help = "Features to include, separated by commas")
    parser.add_argument("-n", "--nucl_seq", action = "store_true", required = False, dest = "nucl_seq", help = "Enable this to include nucleotide sequences of features")
    parser.add_argument("-p", "--prot_seq", action = "store_true", required = False, dest = "prot_seq", help = "Enable this to include protein sequences of CDS")
    return parser.parse_args()  # Returns an instance of the class ArgumentParser

def get_input_filenames(gbks):
    gbk_list = list(gbks)
    if len(gbk_list) == 1 and gbk_list[0].startswith("*"):  # *.gbk
        gbk_list = glob.glob(str(pathlib.Path(".", gbk_list[0])))  # Get names of all GenBank files under the current working directory
    return gbk_list

def prepare_header(args):
    header = ["Contig", "Locus", "Feature", "Start", "End", "Strand", "Pseudo", "Product", "Note"]
    if args.nucl_seq:
        header.append("Nucl_seq")
    if args.prot_seq:
        header.append("Prot_seq")
    return header

def process_records(records, features, args):
    lines = []
    for record in records:
        contig = record.name
        for feature in filter(lambda f: f.type in features, record.features):
            line = process_feature(feature, contig, args)
            lines.append(line)
    return lines

def process_feature(feature, contig, args):
    locus_tag = feature.qualifiers.get("locus_tag", ["unnamed"])[0]
    strand = "+" if feature.strand == 1 else "-"
    is_pseudo = "Y" if "pseudo" in feature.qualifiers or "pseudogene" in feature.qualifiers else "N"
    product = feature.qualifiers.get("product", ["unknown"])[0]
    note = feature.qualifiers.get("note", ["unknown"])[0]

    line = [contig, locus_tag, feature.type, str(feature.location.nofuzzy_start + 1), str(feature.location.nofuzzy_end), strand, is_pseudo, product, note]
    if args.nucl_seq:
        line.append(str(feature.extract(record.seq)))
    if feature.type == "CDS" and args.prot_seq:
        line.append(feature.qualifiers.get("translation", ["unknown"])[0])
    return line

def main():
    args = parse_args()
    gbk_list = get_input_filenames(args.gbks)
    features = args.features.split(",")
    header = prepare_header(args)
    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    for gbk in gbk_list:
        tsv_name = pathlib.Path(args.outdir, pathlib.Path(gbk).stem + ".tsv")
        records = list(SeqIO.parse(gbk, "genbank"))
        lines = process_records(records, features, args)
        
        with open(tsv_name, "w") as tsv:
            tsv.write("\t".join(header) + "\n")
            for line in lines:
                tsv.write("\t".join(line) + "\n")    
                
if __name__ == "__main__":
    main()
