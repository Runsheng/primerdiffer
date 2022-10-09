#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 9/13/16 3:56 PM
# @Author  : Runsheng
# @Email   : Runsheng.lee@gmail.com
# @File    : test_walk_chr.py

# standard library import
import unittest
import os


# self import
from primerdiffer.utils import chr_select,dic2dic,fasta2dic
from primerdiffer.walk_chr import walk_chr_dense, flow_walk_chr
from primerdiffer.logger import log_summary

logger=log_summary()


class TestWalkChr(unittest.TestCase):

    def setUp(self):
        logger.info("****set up****")

        self.sp9_genome = "/t1/ref_BN/cn3_new.fa"
        self.cb4_genome = "/t1/ref_BN/cb5.fa"
        #self.name, self.seq = chr_select(self.cb4_genome, "ChrX", 19424476,19655212)


    def test_walk_chr_dense(self):
        wkdir=os.getcwd()
        pos_str="ChrX:10424476-10655212"
        flow_walk_chr(wkdir=wkdir,
                      genome1=self.cb4_genome,
                      genome2=self.sp9_genome,
                      pos_str=pos_str,
                   interval=10000, jump=2000, out_prefix="primers")

    def tearDown(self):
        logger.info("****tear down****")
        self=None

if __name__ == '__main__':
    unittest.main()