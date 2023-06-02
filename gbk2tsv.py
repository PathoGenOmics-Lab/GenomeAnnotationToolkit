#! /usr/bin/python3
# -*- coding: utf-8 -*-
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
from Bio import SeqIO, SeqFeature
import os
import glob
import sys
import pathlib
from datetime import datetime

__author__ = 'Paula Ruiz-Rodriguez'
__credits__ = ['Paula Ruiz-Rodriguez', 'PathoGenOmics', 'Yu Wan']
__license__ = ['GPL 3', 'Modified from Copyright 2019 Yu Wan (wanyuac@sina.cn)']
__version__ = '0.1.0'
__maintainer__ = 'Paula Ruiz-Rodriguez: @paururo'
__email__ = 'paula.ruiz-rodriguez@uv.es'
__status__ = 'finished'
__lab__ = 'PathoGenOmics (I2SysBio)'

def setup_logging():
    current_date = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')  # Get current date and time
    logger = logging.getLogger(__name__) # Create a custom logger
    logger.setLevel(logging.INFO) # Set the level of this logger
    handler = logging.FileHandler(f'gbk2tsv_{current_date}.log') # Create a file handler
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s') # Create a logging format
    handler.setFormatter(formatter) # Add the formatter to the handler
    logger.addHandler(handler) # Add the handler to the logger  
    return logger

def parse_args():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter,
        description='''
        Transforms GenBank files into tab-separated value files.

        Example commands:
            1. python3 gbk2tsv.py --gbk file.gbk --outdir ./output_directory --features CDS,tRNA,rRNA --nucleotides --protein
            2. python3 gbk2tsv.py --gbk $(ls *.gbk) --outdir . --features CDS,rRNA,tRNA
            3. python3 gbk2tsv.py --gbk file1.gbk file2.gbk --outdir . --features CDS --nucleotides
        
        Requires Biopython
        '''
    )
    parser.add_argument('-g', '--gbk', nargs = '+', type = str, required = True, dest = 'gbks', default = '', help = 'Specify GenBank files as input')
    parser.add_argument('-o', '--outdir', type = str, required = False, dest = 'outdir', default = '.', help = 'Define the directory for output files')
    parser.add_argument('-f', '--features', type = str, required = False, dest = 'features', default = 'CDS,tRNA,rRNA', help = 'Features to include, separated by commas')
    parser.add_argument('-n', '--nucleotides', action = 'store_true', required = False, dest = 'nucleotides', help = 'Enable this to include nucleotide sequences of features')
    parser.add_argument('-p', '--protein', action = 'store_true', required = False, dest = 'protein', help = 'Enable this to include protein sequences of CDS')
    return parser.parse_args()  # Returns an instance of the class ArgumentParser

def get_input_filenames(gbks, logger):
    gbk_list = list(gbks)
    if len(gbk_list) == 1 and gbk_list[0].startswith('*'):  # *.gbk
        gbk_list = glob.glob(str(pathlib.Path('.', gbk_list[0])))  # Get names of all GenBank files under the current working directory
        logger.info(f'Found {len(gbk_list)} GenBank files in the directory')
    return gbk_list

def prepare_header(args):
    header = ['Contig', 'Locus', 'Feature', 'Start', 'End', 'Strand', 'Pseudo', 'Product', 'Note']
    if args.nucleotides:
        header.append('Nucleotide_Seq')
    if args.protein:
        header.append('Protein_Seq')
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
    locus_tag = feature.qualifiers.get('locus_tag', ['unnamed'])[0]
    strand = '+' if feature.strand == 1 else '-'
    is_pseudo = 'True' if 'pseudo' in feature.qualifiers or 'pseudogene' in feature.qualifiers else 'False'
    product = feature.qualifiers.get('product', ['unknown'])[0]
    note = feature.qualifiers.get('note', ['unknown'])[0]
    line = [contig, locus_tag, feature.type, str(feature.location.nofuzzy_start + 1), str(feature.location.nofuzzy_end), strand, is_pseudo, product, note]
    if args.nucleotides:
        line.append(str(feature.extract(record.seq)))
    if feature.type == 'CDS' and args.protein:
        line.append(feature.qualifiers.get('translation', ['unknown'])[0])
    return line

def main():
    logger = setup_logging() # Set up logging
    logger.info('Starting the conversion process')
    
    args = parse_args()
    gbk_list = get_input_filenames(args.gbks, logger)
    
    if args.outdir and not os.path.exists(args.outdir):  
        os.makedirs(args.outdir)
        logger.info(f'Output directory {args.outdir} created')
    
    if (len(gbk_list) == 0):
        logger.error('Invalid --gbk argument: no GenBank file is found.')
        sys.exit('Invalid --gbk argument: no GenBank file is found.')
        
    features = args.features.split(',')
    header = prepare_header(args)
    pathlib.Path(args.outdir).mkdir(parents=True, exist_ok=True)
    
    for gbk in gbk_list:
        logger.info(f'Processing file: {gbk}')
        tsv_name = pathlib.Path(args.outdir, pathlib.Path(gbk).stem + '.tsv')
        records = list(SeqIO.parse(gbk, 'genbank'))
        lines = process_records(records, features, args)
        
        with open(tsv_name, 'w') as tsv:
            tsv.write('\t'.join(header) + '\n')
            for line in lines:
                tsv.write('\t'.join(line) + '\n')
        logger.info(f'Finished processing file: {gbk}')
    
    logger.info('Conversion process completed successfully')   
                
if __name__ == '__main__':
    main()
