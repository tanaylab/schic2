# scHiC 2.0: Sequence and analysis pipeline of single-cell Hi-C datasets #

## Introduction ##

This repository provides the code used to process and analyze single-cell Hi-C libraries in the paper: **Cell-cycle dynamics of chromosomal organisation at single-cell resolution** by Nagano and Lubling et al., Nature, 2017.

It is composed of two main parts (see more details below):

1. **Sequence processing**: relevant code under the _map3c_ folder. It processes a de-multiplexed input fastq paired-end reads files into a list of contacts. 
2. **Single-cell analysis**: relevant code under the _analysis_ folder. It contains the code that builds the data and creates the figures that appear in the paper. 

## Requirements ##
- _Perl_  (with module List::MoreUtils)
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
_misha_ is an R package for genomic analysis, developed by the Tanay lab, and used throughout this repository. To install the package type in your _R_ session:
```
#!r
install.packages("https://schic2.s3.eu-west-1.amazonaws.com

/schic2/misha_3.5.6.tar.gz", repos=NULL)
```

We supply an mm9 genomic database with the genomic and epigenetic data required to run the sequence pipeline and the downstream analysis. Download and unpack [this](https://schic2.s3.eu-west-1.amazonaws.com/schic2_mm9_db.tar.gz) archive (5.5 Gb after unpacking). The _trackdb_ folder in this tar is the root directory of the genomic database, later referred to as _groot_.

## Sequence Processing ##
Processing starts with a pair of fastq files for each cell. The multiplexed fastq files and the cells' indices that allow de-multiplexation are available in GEO under accession [GSE94489](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94489). 
Briefly, Each multiplexed sequencing run comprises 4 FASTQ files: 

1. Forward Sequencing Reads
2. Reverse Sequencing Reads
3. Barcodes 1 (8bps)
4. Barcodes 2 (8bps)

The order of the reads in the FASTQ files correspond to one another. The relevant barcode information can be obtained from the first read in the two barcode files. Reads with unexpected barcodes are written to a single additional file. For more details see the perl script [split_barcodes](https://bitbucket.org/tanaylab/schic2/src/68d7972f64ac2fd32b7c31c5041b39a7176bf14d/map3c/split_barcodes?at=default) (note that it is tailored made to run in Babaraham's cluster environment).


Detailed instructions how to use use the complete pipeline are available in the _map3c_ directory [README](https://github.com/tanaylab/schic2/src/tip/map3c/?at=default)
The final step of the pipeline is to upload the contact map of a cell into the genomic database. We supply all contact maps below to allow you to spare you from rerunning the sequence processing step if you prefer:

|Cells     |Condition  |Batches|Link   |Size |
|:--------:|:---------:|:--------:|:----:|:--------:|
|Haploids |2i |All| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hap_2i_adj_files.tar.gz)|1.2 Gb|
|Haploids |Serum |All| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hap_serum_adj_files.tar.gz)|842 Mb|
|Diploids |2i |1CDU| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hyb_1CDU_adj_files.tar.gz)|461 Mb|
|Diploids |2i |1CDX1| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hyb_1CDX1_adj_files.tar.gz)|511 Mb|
|Diploids |2i |1CDX2| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hyb_1CDX2_adj_files.tar.gz)|618 Mb|
|Diploids |2i |1CDX3| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hyb_1CDX3_adj_files.tar.gz)|618 Mb|
|Diploids |2i |1CDX4| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hyb_1CDX4_adj_files.tar.gz)|779 Mb|
|Diploids |2i |1CDES| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hyb_1CDES_adj_files.tar.gz)|589 Mb|
|Diploids |Serum |1CDS1| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hyb_1CDS1_adj_files.tar.gz)|831 Mb|
|Diploids |Serum |1CDS2| [gz](https://schic2.s3.eu-west-1.amazonaws.com/schic_hyb_1CDS2_adj_files.tar.gz)|1.1 Gb|

Assuming you unpack the archives into a directory pointed by _data_dir_, and you unpack the genomic database into _mm9_db_, run the following code in R to upload contact maps into the genomic database (update _support_sge_ to TRUE if sge is supported):
```
#!r
library(misha)

# define the genomic database we're working on
gdb.init(sprintf("%s/trackdb", mm9_db)) 

# get the list of directories to work on
dirs = list.files(data_dir)

# parse dirs into track names
nms = paste0("scell.nextera.", gsub("-", "_", gsub(".", "_", dirs, fixed=T), fixed=T))

# fends file is the list of GATC fragment ends
fends = sprintf("%s/seq/redb/GATC.fends", mm9_db)

# create tracks directory 
dir.create(sprintf("%s/trackdb/tracks/scell/nextera", mm9_db), showWarnings=F, recursive=T)

# uploading contact files to misha
if (support_sge) {
    # build the commands
    commands = paste0("gtrack.2d.import_contacts(\"", nms, "\", \"\", \"", paste(base_dir, dirs, "adj", sep="/"), "\", fends, allow.duplicates=F)", collapse=",")
    
    # submit jobs to the sge cluster
    res = eval(parse(text=paste("gcluster.run2(",commands,")")))
}
else {
  # upload cell by cell
  for (i in seq_along(nms)) {
    gtrack.2d.import_contacts(nms[i], "", sprintf("%s/%s/adj", data_dir, dirs[i]), fends, allow.duplicates=F)
  }
}
```


#### Translating contact maps into genomic coordinates ####
Each row in the above supplied contact maps represents a single contact. The columns in the files are fend1, fend2 and count, where fend1/2 are the IDs of the interacting fragment ends. You can use the [GATC.fends](http://www.wisdom.weizmann.ac.il/~lubling/schic2/GATC.fends.gz) file (also available in the genomic database supplied above under the seq/redb directory) to translate the fend ids into genomic coordinates. The GATC.fends file maps each fend ID (first column) to chromosome and coordinate (2nd and 3rd columns). 

## Analyse Single-Cell Hi-C datasets ##
Once the contact maps of the cells are uploaded to the genomic database, you can start their analysis. All code and parameter files are supploid under the _analysis_ directory. 

We supply configuration files for the diploids and haploids datasets under _analysis/config_. Make sure to update the paths in those files before you issue any command (look for the "UPDATE REQUIRED" to easily find them).

Run the following commands to build the data and generate all figures (open the R session when you're in the _analysis_ directory):

```
#!r

# load the code
source("paper.r")

# build the tables (need to run once, submits jobs to the sge cluster, takes a while). Supply the diploid (hyb) or haploid (hap) parameters file
sch_build_all(spec_params_fn="config/hyb_2i_params.r")

# After sch_build_all finished, subsequent R sessions can load the tables:
sch_load_tables(spec_params_fn="config/hyb_2i_params.r")

# Once tables are built, finish building the data (need to run once, submits jobs to the sge cluster, takes a while)
build_stats()

# Generate figures
paper_figs()


```

* Metadata information on the cells is available in the features table at [GSE94489](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94489). It contains experimental information on each cell, virtual sorting data and the features that were computed and used in the paper.

## Who do I talk to? ##
For help, please contact yaniv.lubling@weizmann.ac.il
