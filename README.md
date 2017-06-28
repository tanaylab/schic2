# scHiC 2.0: Sequence and analysis pipeline of single-cell Hi-C datasets #

## Introduction ##

This repository provides the code used to process and analyze single-cell Hi-C libraries in the paper: **Cell-cycle dynamics of chromosomal organisation at single-cell resolution** by Nagano and Lubling et al., Nature, 2017.

It is composed of two main parts (see more details below):

1. **Sequence processing**: relevant code under the _map3c_ folder. It processes a de-multiplexed input fastq paired-end reads files into a list of contacts. 
2. **Single-cell analysis**: relevant code under the _analysis_ folder. It contains the code that builds the data and creates the figures that appear in the paper. The _tracdb_ folder in this tar is the root directory of the genomic database, later referred to as _groot_.

## Requirements ##
- _Perl_ 
- [_bowtie2_](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) 
- R with these packages:
    * _misha_ (see below how to install it)
    * [_shaman_](https://bitbucket.org/tanaylab/shaman)
    * _ggplot2_
    * _plyr_
    * _dplyr_
    * _tidyr_
    * _KernSmooth_
    * _RColorBrewer_


## Installing a _misha_ genomic database##
_misha_ is an R package for genomic analysis, developed by the Tanay lab. To install the package type in your _R_ session:
```
#!r
install.packages("http://www.wisdom.weizmann.ac.il/~lubling/schic2/misha_3.5.6.tar.gz", repos=NULL)
```

We supply an mm9 genomic database with the genomic and epigenetic data required to run the sequence pipeline and the downstream analysis. Download and unpack [this](http://www.wisdom.weizmann.ac.il/~lubling/schic2/schic2_mm9_db.tar.gz) archive (5.5 Gb after unpacking). 

## Who do I talk to? ##
For help, please contact yaniv.lubling@weizmann.ac.il