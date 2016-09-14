#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/13/16 2:29 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : primer_check.py


from Bio.Blast.Applications import NcbiblastnCommandline
from cStringIO import StringIO
from Bio.Blast import NCBIXML


def primer_blast(query, db):
    # query = myprimer['PRIMER_LEFT_0_SEQUENCE'] # The sequence
    blastn_cline = NcbiblastnCommandline(db=db, outfmt=5, task="blastn-short") #Blast command
    out, err = blastn_cline(stdin=query)
    blast_records = NCBIXML.read(StringIO(out))  # return is a generator, need a loop to parse the result
    return blast_records


def is_nofalse_primer(blast_records,query,debugmod=False):
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
                    print hsp
                    print hsp.query_end, len(query)
                else:
                    pass
                return False
    # print "pass", hsp
    return True


def primer_check(myprimer, db, primer_number=10, debugmod=False):
    # todo: add login instead of the print function

    '''primer is a return of function primer3.bindings.designPrimers'''
    for i in range(0, primer_number):
        left = myprimer['PRIMER_LEFT_' + str(i) + '_SEQUENCE']
        right = myprimer['PRIMER_RIGHT_' + str(i) + '_SEQUENCE']
        if debugmod:
            print "The %d primer :" % i
            print left, right
        blast_records_l = primer_blast(left,db=db)
        blast_records_r = primer_blast(right,db=db)
        product_size = myprimer['PRIMER_PAIR_' + str(i) + '_PRODUCT_SIZE']

        if is_nofalse_primer(blast_records_l, left, debugmod=debugmod) and is_nofalse_primer(blast_records_r, right,
                                                                                             debugmod=debugmod):
            print "Both pass"
            return (left, right, product_size)
    return 0

