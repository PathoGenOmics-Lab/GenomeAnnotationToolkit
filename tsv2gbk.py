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
        
        for row in tsv_reader:
            if row[0].startswith('#'):
                continue
            feature_type = row[1]
            start = int(row[2]) - 1  # Convert to 0-based index
            end = int(row[3])
            strand = 1 if row[4] == '+' else -1 if row[4] == '-' else 0
            locus_tag = row[5]
            gene = row[6]
            product = row[7]
            db_xrefs = row[8]
            
            qualifiers = {
                'locus_tag': locus_tag,
                'gene': gene,
                'product': product,
                'db_xref': db_xrefs.split(', '),
            }
            if feature_type == 'cds':
                sequence = seq_record.seq[start:end]
                if strand == -1:
                    sequence = sequence.reverse_complement()
                translation = str(sequence.translate(to_stop=False))
                qualifiers['translation'] = translation
                
                gene_qualifiers = {
                    'locus_tag': locus_tag,
                    'gene': gene
                }
                gene_feature = SeqFeature(FeatureLocation(start, end, strand=strand), type='gene', qualifiers=gene_qualifiers)
                seq_record.features.append(gene_feature)
            
            feature = SeqFeature(FeatureLocation(start, end, strand=strand), type=feature_type, qualifiers=qualifiers)
            seq_record.features.append(feature)





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