#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/13/16 2:29 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : primer_check.py


try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

import primer3
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from primerdiffer.general_settings import primer3_general_settings
from primerdiffer.utils import fasta2dic,dic2dic,chr_select, reverse_complement


def my_design_primer(name,seq,primer3_settings=primer3_general_settings):
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



def primer_blast(query, db):
    # query = myprimer['PRIMER_LEFT_0_SEQUENCE'] # The sequence
    blastn_cline = NcbiblastnCommandline(db=db, outfmt=5, task="blastn-short") #Blast command
    out, err = blastn_cline(stdin=query)
    blast_records = NCBIXML.read(StringIO(out))  # return is a generator, need a loop to parse the result
    return blast_records


def filter_hsp(blast_records,query,cutoff_alignlength=16,cutoff_free3=2, debugmod=False):
    """
    filter the hsp, keep only the hsps with align_length <cutoff or free3 <cutoff
    used mainly for insilicon_pcr
    """
    keep=[]
    for n,alignment in enumerate(blast_records.alignments):
        # get all possible alignment position that pass the filter
        if debugmod==True:
            print ("chro is", alignment.hit_def)
        for hsp in alignment.hsps:
            if debugmod==True:
                print ("The subj end is", hsp.sbjct_end)
                print ("The query is", query)
                print (hsp)
                print (hsp.query_end, len(query))
            # get the cutoff
            if hsp.align_length>=cutoff_alignlength or len(query)-hsp.query_end<=cutoff_free3:
                strand=hsp.frame[-1]  # the query is always 1, the target may be -1, 1 is plus and -1 is minus
                keep.append((alignment.hit_def, hsp.sbjct_start, hsp.sbjct_end, strand)) # no end , just one pos
                if debugmod==True:
                    print ("===============Keep==============")

    return keep


def insilicon_pcr(primer_left, primer_right, db, cutoff_alignlength=16, cutoff_free3=2, product_cutoff=2000,
                  debugmod=False):
    """
    para: the left and right primers
    return: a bed-like tuple-list
    """
    possible_product = []

    blast_records_left = primer_blast(primer_left, db)
    blast_records_right = primer_blast(primer_right, db)

    # p_left is a bed like tuple list like [("I", 10000)]
    p_left = filter_hsp(blast_records_left, primer_left, cutoff_alignlength, cutoff_free3)
    p_right = filter_hsp(blast_records_right, primer_right, cutoff_alignlength, cutoff_free3)

    for pl in p_left:  # may need to add a score sys to this function
        for pr in p_right:
            # print pl, pr
            # use only the start to get a approx length, also, the direction for left and right primer should be different
            if pl[0] == pr[0] and abs(pl[1] - pr[1]) <= product_cutoff and pl[-1] * pr[-1] == -1:
                possible_product.append((pl[0], pl[1], pr[1]))

    return possible_product


def product2seqdic(product_l, seqdic):
    """
    use one line from product_l to write the sequence
    """
    outdic={}
    for product_1 in product_l:
        chro, start, end=product_1
        # the product could be forward (as usual) or reverse (F as R and R as F)
        if end>=start:
            name, seq=chr_select(seqdic, chro, start, end)
            outdic[name]=seq
        else:
            name, seq=chr_select(seqdic, chro, end, start)
            outdic[name+"_RC"]=reverse_complement(seq)
    return outdic


def _is_nofalse_primer(blast_records,query,debugmod=False):
    """
    :param blast_records: input a blast record in XML format
    :param query: the query sequence (str)
    :param debugmod: enable the debug mod to see more info
    :return: True or False
    """
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            # so the cutoff is less than 16 match and more than 2 mismatch in 3' region
            if hsp.align_length>=16 or len(query)-hsp.query_end<1:
                if debugmod==True:
                    print (hsp)
                    print (hsp.query_end, len(query))
                else:
                    pass
                return False
    # print "pass", hsp
    return True


def primer_check(myprimer, db1, db2, primer_number=5,
                 cutoff_alignlength=16,cutoff_free3=2, profuct_cutoff=2000,
                 db1_maxhit=1, db2_maxhit=0,
                 debugmod=False):
    '''primer is a return of function primer3.bindings.designPrimers
    db1 and db2 are blastdb,

    CASE1:
    db1 is the fasta genome used to design primer, so the product need to be only 1 (db1_maxhit=1)
    db2 is the fasta genome which should not be amplified, so the product need to be 0 (db2_maxhit=2)

    CASE2:
    db1 is a short sequence used to design primer, products (db1_maxhit) need to be <=1
    db2 is the whole genome, which should has little

    cutoff_alignlength=16,cutoff_free3=2, profuct_cutoff=2000 are used for in silicon PCR
    '''

    for i in range(0, primer_number):
        left = myprimer['PRIMER_LEFT_' + str(i) + '_SEQUENCE']
        right = myprimer['PRIMER_RIGHT_' + str(i) + '_SEQUENCE']
        # designed primer size
        product_size = myprimer['PRIMER_PAIR_' + str(i) + '_PRODUCT_SIZE']

        # the original sequence to detect the false primer
        product_l1=insilicon_pcr(left, right, db1,
                                cutoff_alignlength,cutoff_free3, profuct_cutoff,
                                debugmod=debugmod)
        # the genome used to check false priming
        product_l2=insilicon_pcr(left, right, db2,
                                cutoff_alignlength, cutoff_free3, profuct_cutoff,
                                debugmod=debugmod)
        if debugmod:
            print("The %d primer :" % i)
            print(left, right)
            print("product_l1", product_l1)
            print("produect_l2", product_l2)

        if len(product_l1)<=db1_maxhit and len(product_l2)<=db2_maxhit: # no false primer
            return left, right, product_size # return is a tuple
    return 0


