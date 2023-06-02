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
