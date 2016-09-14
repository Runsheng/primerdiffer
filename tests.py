#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/13/16 3:56 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : tests.py

# standard library import
import unittest


# self import
from utils import chr_select,dic2dic,fasta2dic
from primer_check import primer_check
from walk_chr import my_design_primer, walk_chr_dense

from logger import log_summary
logger=log_summary()


class TestPrimerCheck(unittest.TestCase):

    def setUp(self):
        logger.info("****set up****")
        self.sp9_genome = dic2dic(fasta2dic("./test/sp9pseudo.fa"))
        self.cb4_genome = dic2dic(fasta2dic("./test/cb4.fa"))
        self.name, self.seq = chr_select(self.sp9_genome, "cniII", 700000,740000)

    def tearDown(self):
        logger.info("****tear down****")
        self=None

    def test_1(self):
        myprimer=my_design_primer(name=self.name,seq=self.seq)
        print primer_check(myprimer,db="./test/cb4", debugmod=True)



class TestWalkChr(unittest.TestCase):

    def setUp(self):
        logger.info("****set up****")
        self.sp9_genome = dic2dic(fasta2dic("./test/sp9pseudo.fa"))
        self.cb4_genome = dic2dic(fasta2dic("./test/cb4.fa"))
        self.name, self.seq = chr_select(self.cb4_genome, "X", 19424476,19655212)

        self.seqdic={self.name:self.seq}

    def tearDown(self):
        logger.info("****tear down****")
        self=None

    def test_walk_chr_dense(self):
        walk_chr_dense(self.seqdic, self.name, db="./test/sp9pseudo.fa",
                   interval=4000, jump=400, out_prefix="primers")


if __name__ == '__main__':
    unittest.main()