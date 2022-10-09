#!/usr/bin/env python
#-*- coding: utf-8 -*-
#@Time    : 9/10/2022 5:25 PM
#@Author  : Runsheng
#@File    : ispcr.py

"""
in silico PCR script
"""

import argparse
import os,sys

# self import
from primerdiffer.primer_check import insilicon_pcr, product2seqdic
from primerdiffer.utils import checkblastdb, dic2fasta, dic2dic,fasta2dic


parser=argparse.ArgumentParser()
parser.add_argument("-d", "--wkdir", default=None,
                    help="The dir path contain the file, if not given, use the current dir")
# input
parser.add_argument("-f", "--forward",
                    help="the forward primer sequence")
parser.add_argument("-r", "--reverse",
                    help="the reverse primer sequence")

parser.add_argument("-g", "--genome",
                    help="the fasta file used to design primer")

# primer cutoff
parser.add_argument("--alignlen", type=int, default=16,
                    help="the cutoff of primer min align length as a right hit, default is 16")
parser.add_argument("--free3len",  type=int, default=2,
                    help="the cutoff of primer 3' align length as a right hit, default is 2")
parser.add_argument("--productlen",  type=int, default=2000,
                    help="the cutoff of max product which will be treated as a false priming, default is 2000.")

# output paramenter
parser.add_argument("-o", "--out", default="product.fasta",
                    help="output file, contains all possible amplification regions of this primer pair")

# default handler
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
wkdir=os.getcwd() if args.wkdir is None else args.wkdir

# main
checkblastdb(args.genome) # check blastdb, if no, run makeblastdb

product_l=insilicon_pcr(primer_left=args.forward,
              primer_right=args.reverse,
              db=args.genome,
              cutoff_alignlength=args.alignlen,
              cutoff_free3=args.free3len,
              product_cutoff=args.productlen
              )
seqdic=dic2dic(fasta2dic(args.genome))
outdic=product2seqdic(product_l,seqdic)
dic2fasta(outdic, args.out)





