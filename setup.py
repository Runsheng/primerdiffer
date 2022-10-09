from setuptools import setup
from primerdiffer import __version__

setup(
    name='primerdiffer',
    version=__version__,
    packages=['', "primerdiffer", "bin"],
    url='https://github.com/Runsheng/primerdiffer',
    license='GPL-2',
    author='runsheng',
    author_email='runsheng.lee@gmail.com',
    description='batch design of primers for genotyping',
    install_requires = ["primer3-py>=0.6.1",
                        "biopython>=1.78"],
    scripts = ['bin/primerdesign.py',
               'bin/ispcr.py'
                ]
)

