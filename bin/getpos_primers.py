#!/usr/bin/env python
#-*- coding: utf-8 -*-
#@Time    : 17/7/2023 10:18 pm
#@Author  : Runsheng
#@File    : getpos_primers.py

"""
get the position using ispcr function for all output from primerdiffer.py
run example:
cd /t1/ref_BN
getpos_primers.py -f primer10.txt -g cb5.fa -o primer10_pos.txt
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
parser.add_argument("-f", "--file",
                    help="the 4 col table generated by primerdesign.py")

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
parser.add_argument("-o", "--out", default="primerpos.txt",
                    help="output file, contains possible amplification regions of this primer pair")

# default handler
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
wkdir=os.getcwd() if args.wkdir is None else args.wkdir

# main
checkblastdb(args.genome) # check blastdb, if no, run makeblastdb
seqdic=dic2dic(fasta2dic(args.genome))

def read_primerdiffer_out(filename):
    lines=[]
    with open(filename, "r") as f:
        for line in f.readlines():
            line_l=line.strip().split("\t")
            lines.append(line_l)
    return lines

lines=read_primerdiffer_out(args.file)

lines_new=[]
for line_l in lines:
    name_raw, forward, reverse, length_str=line_l

    product_l=insilicon_pcr(primer_left=forward,
              primer_right=reverse,
              db=args.genome,
              cutoff_alignlength=args.alignlen,
              cutoff_free3=args.free3len,
              product_cutoff=args.productlen
              )

    outdic=product2seqdic(product_l,seqdic)
    name_pos=list(outdic.keys())[0]
    seq=outdic[name_pos]
    line_one=[name_raw, forward, reverse, length_str, name_pos, seq]
    lines_new.append(line_one)

with open(args.out, "w") as fw:
    for line_one in lines_new:
        fw.write("\t".join(line_one))
        fw.write("\n")