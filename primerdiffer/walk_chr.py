#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/13/16 3:23 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : walk_chr.py

# third part import
import primer3

# self import
from utils import chr_select
from primer_check import primer_check
from general_settings import primer3_general_settings


def my_design_primer(name,seq,primer3_settings=primer3_general_settings):
    # todo: to make the setting for the primer number easier
    """
    general wrapper for primer3-py
    :param name: name for the sequence
    :param seq: string for sequence in UPPER case
    :param primer3_settings: general setting for primer_design
    :return: the dict storing the primer pairs specific for seq
    """
    seq_args = {'SEQUENCE_ID': name,
                'SEQUENCE_TEMPLATE':seq}
    myprimer=primer3.bindings.designPrimers(seq_args,primer3_settings)
    return myprimer


def walk_chr_dense(genome, chro, db,
                   interval=500000, jump=4000, out_prefix="primers"):
    """
    :param genome: genome is a dict in name:seq,
    :param chro:
    :param interval: interval is the region begins to pick primers
    :param jump: jump is the region to jump inside intervals, the smaller, more sites to check
    :param out_prefix:
    :return:
    """
    primer_dict = {}
    f_out = open(out_prefix + "_" + chro + ".txt", "w")

    n = len(genome[chro]) / interval
    i = 0
    while i < n:
        offset = 0
        while offset <= (interval / 2) / jump:
            name, seq = chr_select(genome, chro, i * interval + offset * jump,
                                   i * interval + (offset + 1) * jump)
            # print name,seq
            if "N" in seq.upper():
                offset += 1
            else:
                myprimer = my_design_primer(name=name, seq=seq)
                primer_status = primer_check(myprimer,db)
                if primer_status == 0:
                    offset += 1
                else:
                    offset = interval / jump  # just >(interval/2)/jump, used to indicate the sucess of primer finding
                    left, right, product_size = primer_status
        i += 1

        if offset == interval / jump:
            print("Get primer in %s!" % name)
            primer_dict[name] = (left, right, product_size)
            f_out.write(name + "\t" + left + "\t" + right + "\t" + str(product_size) + "\n")
        if offset == (interval / 2) / jump + 1:
            print ("Can not get primer in %s!" % name)

    f_out.close()
    return primer_dict
