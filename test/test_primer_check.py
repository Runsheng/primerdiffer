#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/14/16 4:59 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : test_primer_check.py


import unittest

# self import
from primerdiffer.utils import chr_select,dic2dic,fasta2dic
from primerdiffer.primer_check import *


class TestPrimerCheck(unittest.TestCase):

    def setUp(self):
        self.sp9_genome = "/t1/ref_BN/cn3_new.fa"
        self.cb4_genome = "/t1/ref_BN/cb5.fa"
        self.cb4=dic2dic(fasta2dic(self.cb4_genome))

    def test_design_primer(self):
        name, seq = chr_select(self.cb4, 'ChrI', 231537 - 1000, 231551 + 1000)
        myprimer = my_design_primer(name=name, seq=seq)
        for k,v in myprimer.items():
            print (k,v)

    def test_primer_blast(self):
        blast_records = primer_blast("gggtgagaatagagtggtgg",
                                     db=self.cb4_genome)
        for i in blast_records.alignments:
            print(i)

    def test_filter_hsp(self):
        blast_records = primer_blast('gggtgagaatagagtggtgg',
                                     db=self.cb4_genome)
        keep = filter_hsp(blast_records, 'gggtgagaatagagtggtgg', debugmod=True)
        print(keep)


    def test_insilicon_pcr(self):
        aa = insilicon_pcr(primer_left='gcactttcatgtccctcaac',
                           primer_right='cactctattctcaccccacc',
                           db=self.cb4_genome)
        ab = insilicon_pcr("gttgagggacatgaaagtgc", 'ctaggcgaaagaggtcacat',
                           db=self.cb4_genome)

        cc= insilicon_pcr("gttgagggacatgaaagtgc", 'aaaaaaagctctctctggg',
                           db=self.cb4_genome)
        dd=insilicon_pcr("AGATTGACCGAATTCAGCCT", "TCCGTTTTGTAGAGCTCCTC",
                         db=self.sp9_genome)
        d2=insilicon_pcr("AGATTGACCGAATTCAGCCT", "TCCGTTTTGTAGAGCTCCTC",
                         db=self.cb4_genome)
        print(aa,ab,cc,dd, d2)

    def test_primer_check(self):
        name, seq = chr_select(self.cb4, 'ChrI', 231537 - 1000, 231551 + 1000)
        myprimer = my_design_primer(name=name, seq=seq)

        # design the primer using
        print( "GET PRIMER:", primer_check(myprimer,
                     db1=self.sp9_genome,
                     db2=self.cb4_genome,
                     debugmod=True)
               )

        print( "GET PRIMER:", primer_check(myprimer,
                            db1=self.cb4_genome,
                            db2=self.sp9_genome,
                            debugmod=True)
               )

    def test_product2seqdic(self):
        aa = insilicon_pcr(primer_left='gcactttcatgtccctcaac',
                           primer_right='cactctattctcaccccacc',
                           db=self.cb4_genome)
        print(product2seqdic(aa, self.cb4))


    def tearDown(self):
        self=None

if __name__ == '__main__':
    unittest.main()