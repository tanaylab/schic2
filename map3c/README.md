## Running the sequence processing pipeline ##

### Before you start ###
Define the following environment variables:

1. _PIPELINE_HOME_ : should point to your repository directory (parent of _map3c_ directory).
2. _PERL5LIB_      : should include the repository directory.
3. _MM9_DB_ : path to the [mm9 genomic database](https://schic2.s3.eu-west-1.amazonaws.com/schic2/schic2_mm9_db.tar.gz) directory you downloaded and unpacked.

We assume you downloaded the FASTQ files from GEO, accession [GSE94489](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE94489) and de-multiplexed them such that each cell has a pair of sequence files (paired-end fastq files) with the cell's indices in their names. You should de-multiplex each fastq into a separate directory to prevent name clashing between cells from different batches but with the same indices.

### Step 1: Generate mappability track ###
To generate the mappability track we chop the genome into short sequences, map those sequences back to the genome and count the fraction of mapped reads in each binsize.

**Note:** The genomic database you've downloaded already contains the mapability track, which was produced with the parameter file _config/mm9_mapa.conf_.

If you wish to produce another mappability track, first update the paths to your bowtie2 in _config/mm9_mapa.conf_., then launce an R session from within the _map3c_ directory and type:

```
#!r
source("TG3C/imp3c.r")
library(misha)
gtrack.create_mapab_track(params_fn="config/mm9_mapa.conf")
```

### Step 2: Generate the restriction enzyme fragments tracks ###
To generate the restriction enzyme tables and tracks (redb) we scan the genome with the restriction sequence and generate tables of fragments and the matching fragment ends.

**Note:** The genomic database you've downloaded already contains the redb data, which was produced with the parameter file _config/mm9_redb.conf_.

To reproduce the redb, launce an R session from within the _map3c_ directory and type:

```
#!r
source("TG3C/imp3c.r")

# re_seq is the recognition site, GATC for MboI or DpnII enzymes.
gtrack.create_redb_tracks(re_seq="GATC", params_fn="config/mm9_redb.conf")
```

### Step 3: Processing single-cell dataset
Processing of each cell starts with the de-multiplexed fastq file and ends with computing the contact map and uploading it to the genomic database as a 2D tack. The processing is managed by configuration files. Each configuration file can be submitted as a job to allow parallel processing.
Here's a sample config file with annotations:

```

# cell names to process in the current job. 
TG3C.3C_indices=326,327,328

# For each TG3C.3C_indices specified above add the directory where its fastq filesare and the regular expression to match its fastq in that directory (use the pair of indices unique for that cell)
TG3C.3C_dir_326=/net/mraid14/export/tgdata/db/tgdb/mm9/rawdata/scell_hic/cells_hyb_apr_2016/Sample_3222/
TG3C.3C_fn_regexp_326=GGACTCCT_GTAAGGAG
TG3C.3C_dir_327=/net/mraid14/export/tgdata/db/tgdb/mm9/rawdata/scell_hic/cells_hyb_apr_2016/Sample_3222/
TG3C.3C_fn_regexp_327=GGACTCCT_TATCCTCT
TG3C.3C_dir_328=/net/mraid14/export/tgdata/db/tgdb/mm9/rawdata/scell_hic/cells_hyb_apr_2016/Sample_3222/
TG3C.3C_fn_regexp_328=GGACTCCT_TCTCTCCG

# shared definitions
#include ${PIPELINE_HOME}/map3c/config/scell_shared.conf

# cell names prefix. Should be either NXT (for haploid 2i cells), NST (for haploid serum), 1CDU/1CDX1/1CDX2/1CDX3/1CDX4/1CDES (for diploid 2i) or 1CDS1/1CDS2 (for dipliod serum)
TG3C.3C_exp_nm=1CDX3

```

Use the script _TG3C/split_pairs_to_cfgs.pl_ to generate these configuration files. Run it for each sample separetely (to avoid index clashing). It read from the standard input lines containing the unique regexp of a cell and the cell number (tab delimited) and generate configuration files according to the supplied parameters. 

Example run (replace capital parameters with real parameters):

```
cat CELLS_REGEXP_NAMES_PAIRS.TXT | perl ${PIPELINE_HOME}/map3c/TG3C/split_pairs_to_cfgs.pl INPUT_FASTQ_DIR OUTPUT_FILE_PREFIX 50 1CDX1
```

Update the paths to the bowtie2 executable and the mm9 indices in the map3c/config/scell_shared.conf file.

Once all config files are generated, you can issue a single command per config file. From within R (opened under the _map3c_ directory), type (replace CONFIG_FILE with the actual value):

```
#!r
source("TG3C/imp3c.r")
gtrack.2d.create_from_3Cseq(newtrack_name="scell.nextera", track_desc="schic track", groot=sprintf("%s/trackdb", Sys.getenv("MM9_DB")), params_fn=CONFIG_FILE)

```

To issue this command from the shell (and allow you to submit it as a job) you can use the _TG3C/call.r_ script we supply:
```
#!bash
${PIPELINE_HOME}/map3c/TG3C/call.r $PIPELINE_HOME/map3c TG3C/imp3c.r gtrack.2d.create_from_3Cseq newtrack_name=scell.nextera track_desc=schic groot=${MM9_DB}/trackdb/ params_fn=CONFIG_FILE


```
