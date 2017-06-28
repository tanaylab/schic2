# scHiC 2.0: Sequence and analysis pipeline of single-cell Hi-C datasets #

## Introduction ##

This repository provides the code used to process and analyze single-cell Hi-C libraries in the paper: **Cell-cycle dynamics of chromosomal organisation at single-cell resolution** by Nagano and Lubling et al., Nature, 2017.


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


#### Installing misha package:
Type in your _R_ session:
```
#!r
install.packages("http://www.wisdom.weizmann.ac.il/~lubling/schic2/misha_3.5.6.tar.gz", repos=NULL) # Download and install misha package
```

## Who do I talk to? ##
For help, please contact yaniv.lubling@weizmann.ac.il