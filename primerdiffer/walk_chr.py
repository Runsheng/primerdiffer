#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/13/16 3:23 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : walk_chr.py

# third part import

# self import
import os

from primerdiffer.primer_check import my_design_primer, primer_check
from primerdiffer.utils import dic2dic, fasta2dic, chr_select, tuple_to_pos_str, pos_str_to_tuple, checkblastdb


def walk_chr_dense(genome, chro, start, end, db1, db2,
                   cutoff_alignlength=16, cutoff_free3=2, product_cutoff=2000,
                   db1_maxhit=1, db2_maxhit=0,
                   interval=500000, jump=4000, out_prefix="primers", debugmod=False):
    """
    :param genome: genome is a dict in name:seq,
    :param chro:
    :param interval: interval is the region begins to pick primers
    :param jump: jump is the region to jump inside intervals, the smaller, more sites to check
    :param out_prefix:
    :return:
    """
    primer_dict = {}
    pos_str= tuple_to_pos_str((chro, start, end))
    f_out = open(out_prefix + "_" + pos_str + ".txt", "w")

    n = len(genome[chro][start:end]) / interval
    i = 0
    while i < n:
        offset = 0
        while offset <= (interval / 2) / jump:
            name, seq = chr_select(genome, chro, start+i * interval + offset * jump,
                                   start+i * interval + (offset + 1) * jump)
            print ("Design primers for ", name)
            if "N" in seq.upper():
                offset += 1
            else:
                myprimer = my_design_primer(name=name, seq=seq)
                primer_used = primer_check(myprimer,db1,db2,
                                           cutoff_alignlength, cutoff_free3, product_cutoff,
                                           db1_maxhit, db2_maxhit,
                                           debugmod=debugmod)
                if primer_used == 0:
                    offset += 1
                else:
                    offset = interval / jump  # just >(interval/2)/jump, used to indicate the sucess of primer finding
                    left, right, product_size = primer_used
        i += 1

        if offset == interval / jump:
            print("Get primer in %s!" % name)
            primer_dict[name] = (left, right, product_size)
            f_out.write(name + "\t" + left + "\t" + right + "\t" + str(product_size) + "\n")
        if offset == (interval / 2) / jump + 1:
            print ("Can not get primer in %s!" % name)

    f_out.close()
    return primer_dict


def flow_walk_chr(wkdir, genome1, genome2, pos_str,
                  cutoff_alignlength=16, cutoff_free3=2, product_cutoff=2000,
                  db1_maxhit=1, db2_maxhit=0,
                  interval=4000, jump=400, out_prefix="primers"):
    """

    genome1: the genome fasta file used to design primers
    genome2: the genome fasta used to blast against, the primer should not amplify any product in genome2

    pos_str: chro:start-end, same as the input for igv or ucsc web browser

    """
    # can add wkdir is os.getcwd() later
    os.chdir(wkdir)

    # first to test if the blastdb for genome2 is there, if not, create one
    checkblastdb(genome1)
    checkblastdb(genome2)

    # main
    g1=dic2dic(fasta2dic(genome1))
    chro, start, end = pos_str_to_tuple(pos_str)

    primer_dict=walk_chr_dense(genome=g1, chro=chro, start=start, end=end,
                               db1=genome1,db2=genome2,
                               cutoff_alignlength=cutoff_alignlength,
                               cutoff_free3=cutoff_free3,
                               product_cutoff=product_cutoff,
                               db1_maxhit=db1_maxhit,
                               db2_maxhit=db2_maxhit,
                            interval=interval, jump=jump, out_prefix=out_prefix)
    print(primer_dict)
