#!/usr/bin/env python
#-*- coding: utf-8 -*-
#@Time    : 9/10/2022 4:30 PM
#@Author  : Runsheng
#@File    : primerdesign.py

"""
The main script used to run cmd orders
"""

import argparse
import os
import sys

#currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
#parentdir = os.path.dirname(os.path.dirname(currentdir))
#sys.path.insert(0,parentdir)

# self import
from primerdiffer.walk_chr import flow_walk_chr, primer3_general_settings

parser=argparse.ArgumentParser()
parser.add_argument("-d", "--wkdir", default=None,
                    help="The dir path contain the file, if not given, use the current dir")
# input
parser.add_argument("-g1", "--genome1",
                    help="the fasta file used to design primer")
parser.add_argument("-g2", "--genome2",
                    help="the fasta file used to check false priming")
parser.add_argument("-pos", "--position",default=None,
                    help="position on genome1 to design primer, use string with IGV format, like ChrX:1956230-1976220")

# primer cutoff
parser.add_argument("--alignlen", type=int, default=16,
                    help="the cutoff of primer min align length as a right hit, default is 16")
parser.add_argument("--free3len",  type=int, default=2,
                    help="the cutoff of primer 3' align length as a right hit, default is 2")
parser.add_argument("--productlen",  type=int, default=2000,
                    help="the cutoff of max product which will be treated as a false priming, default is 2000.")

# hit cutoff
parser.add_argument("-h1","--hit1", type=int, default=1,
                    help="the cutoff of max number of in-silicon PCR product which can be found in genome1. default is 1")
parser.add_argument("-h2","--hit2",  type=int, default=0,
                    help="the cutoff of max number of in-silicon PCR product which can be found in genome2, default is 0")

# run parameter, how dense the primer should be
parser.add_argument("-i","--interval", type=int, default=5000,
                    help="interval is the region begins to pick primers, default is 5000. If 5k is the unit, will pick one primer each 5k")
parser.add_argument("-j","--jump",  type=int, default=500,
                    help="jump is the region to jump inside intervals if the prvious interval can not generate a valid primer"
                         ", the smaller, more sites to check. Default is 500.")
# output paramenter
parser.add_argument("--prefix", default="primers",
                    help="prefix of output file, default is primers")

# add user settings
parser.add_argument("--primer3config", default=None,
                    help="the config file for the primer3 ")

# default handler
args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
wkdir=os.getcwd() if args.wkdir is None else args.wkdir


# deal with config
primer3_setting = primer3_general_settings
if args.primer3config is None:
    pass
else: # update the primer3 setting from the given parameters
    with open(args.primer3config, "r") as f:
        primer3_user_setting = eval(f.read()) # read the config as a dict
        primer3_setting.update(primer3_user_setting)


flow_walk_chr(wkdir=wkdir,
              genome1=args.genome1,
              genome2=args.genome2,
              pos_str=args.position,
              cutoff_alignlength=args.alignlen,
              cutoff_free3=args.free3len,
              product_cutoff=args.productlen,
              db1_maxhit=args.hit1,
              db2_maxhit=args.hit2,
              interval=args.interval,
              jump=args.jump,
              out_prefix=args.prefix,
              primer3_settings=primer3_setting)
