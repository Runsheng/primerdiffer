#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/13/16 3:56 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : test_walk_chr.py

# standard library import
import unittest


# self import
from primerdiffer.utils import chr_select,dic2dic,fasta2dic
from primerdiffer.walk_chr import walk_chr_dense

from primerdiffer.logger import log_summary
logger=log_summary()


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