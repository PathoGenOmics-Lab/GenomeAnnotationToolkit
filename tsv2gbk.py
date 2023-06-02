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



def parse_args():
    parser = argparse.ArgumentParser(description='Convert a TSV file to a GenBank file')
    parser.add_argument('-tsv', dest='tsv_file', required=True, help='Input TSV file')
    parser.add_argument('-fasta', dest='fasta_file', required=True, help='Input FASTA file')
    parser.add_argument('-o', dest='gbk_file', required=True, help='Output GenBank file')
    parser.add_argument('-molecule_type', dest='molecule_type', default='DNA', help='Molecule type (default: DNA)')
    
    return parser.parse_args()  # Returns an instance of the class ArgumentParser