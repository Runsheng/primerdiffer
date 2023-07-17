from setuptools import setup
from primerdiffer import __version__

with open('readme.md', "r") as f:
    LONG_DESC = f.read()

setup(
    name='primerdiffer',
    version=__version__,
    packages=['', "primerdiffer", "bin"],
    url='https://github.com/Runsheng/primerdiffer',
    license='GPL-2',
    author='runsheng',
    author_email='runsheng.lee@gmail.com',
    description='batch design of primers for genotyping two or more haplotypes in hybrids',
    long_description=LONG_DESC,
    long_description_content_type='text/markdown',
    install_requires = ["primer3-py>=0.6.1",
                        "biopython>=1.78"],
    scripts = ['bin/primerdesign.py',
               'bin/ispcr.py',
               'bin/getpos_primers.py'
                ]
)

