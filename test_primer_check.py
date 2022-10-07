#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/14/16 4:59 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : test_primer_check.py


import unittest

# self import
from utils import chr_select,dic2dic,fasta2dic
from primer_check import primer_check
from walk_chr import my_design_primer, walk_chr_dense

# set logger
from logger import log_summary
logger=log_summary()


class TestPrimerCheck(unittest.TestCase):

    def setUp(self):
        logger.info("****set up %s ****" % self.__name__)
        self.sp9_genome = dic2dic(fasta2dic("./test/sp9pseudo.fa"))
        self.cb4_genome = dic2dic(fasta2dic("./test/cb4.fa"))
        self.name, self.seq = chr_select(self.sp9_genome, "cniII", 700000,740000)


    def test_1(self):
        myprimer=my_design_primer(name=self.name,seq=self.seq)
        print (primer_check(myprimer,db="./test/cb4", debugmod=True))

    def tearDown(self):
        logger.info("****tear down %s ****" % self.__str__)
        self=None

if __name__ == '__main__':
    unittest.main()