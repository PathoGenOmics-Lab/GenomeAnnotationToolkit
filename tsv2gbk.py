#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import csv
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Data import CodonTable

def tsv_to_gbk(tsv_file, fasta_file, gbk_file, molecule_type='DNA'):
    with open(tsv_file, 'r') as tsv:
        tsv_reader = csv.reader(tsv, delimiter='\t')
        
        for record in SeqIO.parse(fasta_file, 'fasta'):
            fasta_name = fasta_file.split('/')[-1].split('.')[0]
            seq_record = SeqRecord(Seq(str(record.seq)), id=fasta_name, description=record.description)

        seq_record.annotations['molecule_type'] = molecule_type



def parse_args():
    parser = argparse.ArgumentParser(description='Convert a TSV file to a GenBank file')
    parser.add_argument('-tsv', dest='tsv_file', required=True, help='Input TSV file')
    parser.add_argument('-fasta', dest='fasta_file', required=True, help='Input FASTA file')
    parser.add_argument('-o', dest='gbk_file', required=True, help='Output GenBank file')
    parser.add_argument('-molecule_type', dest='molecule_type', default='DNA', help='Molecule type (default: DNA)')
    
    return parser.parse_args()  # Returns an instance of the class ArgumentParser

def main():
    args = parse_args()

if __name__ == '__main__':
    main()