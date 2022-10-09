# primerdiffer: batch design primers for chromosomal level genotyping 
![PyPI](https://img.shields.io/pypi/v/primerdiffer?color=green)

## Installation:
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
- Greedy design primers for a region in genome1 and make a specifiy check using genome2
- The dis-similarity between genome1 and genome2 >5%. 
```bash
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

Use C. nigoni and C.briggsae genome as example:
```bash
# the C. nigoni genome is cn3_new.fa and C. briggsae genome is cb5.fa
# design C. briggsae unique primer, which would not amplify any region in C. nigoni, 
# in the region of ChrX:12881200-15106660
# inside every 4kb interval
primerdesign.py -g1 cb5.fa -g2 cn3_new.fa -pos "ChrX:12881200-15106660" --interval 4000
head primers_ChrX\:12881200-15106660.txt
#ChrX:12881200-12881700	GATCCAAAACATGAGTGGCC	CGAGATCATTGGCTCAAAGT	287
#ChrX:12886200-12886700	GTTTTCTCTTCAAGTGCCCG	CTCCCACATCTTGTAGGTCC	416
#ChrX:12891200-12891700	GTAGATGCTGTTGAGGCTCT	CGAGTGGGACATTGTCAGTA	299
#ChrX:12896200-12896700	GGCGCATTATACGAAGCTTT	TTCCTGCTGCCAGATAGAAG	362
```


Use in silico PCR to get the position and the product of the primer
```bash
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


    