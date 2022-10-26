# primerdiffer: batch design primers for chromosomal level genotyping 
![PyPI](https://img.shields.io/pypi/v/primerdiffer?color=green)

## Installation:
The package worked with python version >=3.4.
Only tested in linux x64 system.

python package:
- primer3-py>=0.6.1
- biopython>=1.7.8

other program in your $PATH:
- ncbi-blast

example code to install the packages under python3 with pip and conda
```bash
pip install primerdiffer # will also install primer3-py and biopython  
conda install -c bioconda blast # install ncbi blast, which is not included in pip installation
```

## Case example: primerdesign.py for primerdiffer
Design genome-wide specific primers for two species/sub-species/divergent sequences:
- Greedy design primers for a region in genome1 and make a specificity check using genome2
- The dis-similarity between genome1 and genome2 >= 5%. 
```
usage: primerdesign.py [-h] [-d WKDIR] [-g1 GENOME1] [-g2 GENOME2] [-pos POSITION] [--alignlen ALIGNLEN]
                       [--free3len FREE3LEN] [--productlen PRODUCTLEN] [-h1 HIT1] [-h2 HIT2]
                       [-i INTERVAL] [-j JUMP] [--prefix PREFIX]

optional arguments:
  -h, --help            show this help message and exit
  -d WKDIR, --wkdir WKDIR
                        The dir path contain the file, if not given, use the current dir
  -g1 GENOME1, --genome1 GENOME1
                        the fasta file used to design primer
  -g2 GENOME2, --genome2 GENOME2
                        the fasta file used to check false priming
  -pos POSITION, --position POSITION
                        position on genome1 to design primer, use string with IGV format, like
                        ChrX:1956230-1976220
  --alignlen ALIGNLEN   the cutoff of primer min align length as a right hit, default is 16
  --free3len FREE3LEN   the cutoff of primer 3' align length as a right hit, default is 2
  --productlen PRODUCTLEN
                        the cutoff of max product which will be treated as a false priming, default is
                        2000.
  -h1 HIT1, --hit1 HIT1
                        the cutoff of max number of in-silicon PCR product which can be found in
                        genome1. default is 1
  -h2 HIT2, --hit2 HIT2
                        the cutoff of max number of in-silicon PCR product which can be found in
                        genome2, default is 0
  -i INTERVAL, --interval INTERVAL
                        interval is the region begins to pick primers, default is 5000. If 5k is the
                        unit, will pick one primer each 5k
  -j JUMP, --jump JUMP  jump is the region to jump inside intervals if the prvious interval can not
                        generate a valid primer, the smaller, more sites to check. Default is 500.
  --prefix PREFIX       prefix of output file, default is primers
```

Use _C. nigoni_ and _C.briggsae_ genomes as example. The two fasta files can be downloaded separately 
from [cb4.fa](https://github.com/Runsheng/cbgenome/releases/download/cb5pre_cn3pre/cb5.fa.gz) and 
[cn3_new.fa](https://github.com/Runsheng/cbgenome/releases/download/cb5pre_cn3pre/cn3_new.fa.gz). 

The _C. nigoni_ genome is cn3_new.fa and _C. briggsae_ genome is cb5.fa. To design _C. briggsae_ unique primer, 
which would not amplify any region in _C. nigoni_, and amplify only one region in _C. briggsae_. 
The targeted region for C. briggsae is ChrX:12881200-15106660 (-pos),
one primer is designed for every 4kb interval (--interval).
```
primerdesign.py -g1 cb5.fa -g2 cn3_new.fa -pos "ChrX:12881200-15106660" --interval 4000

# check the result in file "primers_ChrX:12881200-15106660.txt"
head primers_ChrX\:12881200-15106660.txt
#ChrX:12881200-12881700	GATCCAAAACATGAGTGGCC	CGAGATCATTGGCTCAAAGT	287
#ChrX:12886200-12886700	GTTTTCTCTTCAAGTGCCCG	CTCCCACATCTTGTAGGTCC	416
#ChrX:12891200-12891700	GTAGATGCTGTTGAGGCTCT	CGAGTGGGACATTGTCAGTA	299
#ChrX:12896200-12896700	GGCGCATTATACGAAGCTTT	TTCCTGCTGCCAGATAGAAG	362
```


Use in silico PCR to get the position and the product of the primer
```
usage: ispcr.py [-h] [-d WKDIR] [-f FORWARD] [-r REVERSE] [-g GENOME] [--alignlen ALIGNLEN]
                [--free3len FREE3LEN] [--productlen PRODUCTLEN] [-o OUT]

optional arguments:
  -h, --help            show this help message and exit
  -d WKDIR, --wkdir WKDIR
                        The dir path contain the file, if not given, use the current dir
  -f FORWARD, --forward FORWARD
                        the forward primer sequence
  -r REVERSE, --reverse REVERSE
                        the reverse primer sequence
  -g GENOME, --genome GENOME
                        the fasta file used to design primer
  --alignlen ALIGNLEN   the cutoff of primer min align length as a right hit, default is 16
  --free3len FREE3LEN   the cutoff of primer 3' align length as a right hit, default is 2
  --productlen PRODUCTLEN
                        the cutoff of max product which will be treated as a false priming, default is
                        2000.
  -o OUT, --out OUT     output file, contains all possible amplification regions of this primer pair
```

## Case example: in silico PCR product

```bash
ispcr.py -f gcactttcatgtccctcaac -r cactctattctcaccccacc -g cb5.fa -o ispcr.fa
head ispcr.fa
#>ChrI:230699-231076_RC
#GCACTTTCATGTCCCTCAACCAGTCGTTTTTCCTTACCTCTCCCCTTCCTTTTTTCCCCCTCCCAGATGACGTCACCCATCTGTCC
#ACCCTTCTAACGGTCCCCCCCACCATCTGCATGGTGTCCTCGGGGGTGAACAGCTGCACATTTATTGTTCCCTTCTATTCCCCCCT
#CCTCGGTCATCGCGGTTTTATCCCCGCCGTCCATTTGACCATTCTTTGCGTCTCTTCCCTCTCTCTCTCTCTCTCTCTCTCTCTCT
#CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTTCTTTTTCTTCTAGGCGAAAGAGGTCACATGGAAGAGAAGAGGATG
#ATGATGATGATGATGGTGGGGTGAGAATAGAGT
```

## Roadmap for other functions:
1. To use user-provided primer parameters.Primer design parameter now is fine-tuned for general purpose PCR, which can be found in "general_settings.py".This file may need be modified to generate primers for specific purpose PCR like real-time qPCR.
2. To update the RFLP method for primer design to differ sequences with almost identical sequence.
3. To update the primer design using VCF file.

    