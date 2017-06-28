source("analyzeScHiC.r")

options(error=recover)


#############
wrapper_sch_calc_tad_exact_borders_insu_ds_trian <- function(chroms = paste0("chr", sch_chroms), n_ds=8, rebuild=F, scale=2e5, res=0, discard_below=1000, name="_tad_borders", ins_track=ins_tn, filter_iter=0)
{
  b = .sch_get_pool_tad_borders(scale=scale, min_diag=discard_below, ins_track=ins_track, filter_iter=filter_iter)
  b_s = split(b$start, as.character(b$chrom))

  commands = paste("sch_calc_insu_ds_trian(\"", chroms, "\", coords=b_s[[\"", chroms, sprintf("\"]], coords_nm=\"%s\", n_ds=%d, rebuild=%s, scale=%d, res=%d, return_data=F, discard_below=%f)", name, n_ds, rebuild, scale, res, discard_below), sep="", collapse=",")
  print(commands)
  res=eval(parse(text=paste("gcluster.run(",commands,", R=\"/net/mraid14/export/data/users/eladch/tools/CO6/R/3.2.1/bin/R\", opt.flags=\"-q all.q@@dell6220-128g\")")))
  res
}

######
wrapper_sch_calc_tad_exact_borders_insu_ds_trian_discard_unear <- function(chroms=paste0("chr", sch_chroms), rebuild=F)
{
  wrapper_sch_calc_tad_exact_borders_insu_ds_trian(chroms=chroms, n_ds=8, rebuild=rebuild, scale=3e5, res=0, discard_below=sch_remove_near_dist, name=sprintf("_tad_borders_gt%s", n2str(sch_remove_near_dist)))
}



####
wrapper_sch_calc_dixon_borders_insu_ds_trian <- function(chroms = paste0("chr", sch_chroms), n_ds=8, rebuild=F, scale=3e5, res=0, discard_below=sch_remove_near_dist, name=sprintf("_Dixon_tad_borders_gt%s", n2str(sch_remove_near_dist)), max_border_size=80e3)
{
  b = read.delim(dixon_borders_fn, header=T)
  b = b[b$end - b$start <= max_border_size, ]
  b = intervals.centers(b)
  

  b_s = split(b$start, as.character(b$chrom))

  commands = paste("sch_calc_insu_ds_trian(\"", chroms, "\", coords=b_s[[\"", chroms, sprintf("\"]], coords_nm=\"%s\", n_ds=%d, rebuild=%s, scale=%d, res=%d, return_data=F, discard_below=%f)", name, n_ds, rebuild, scale, res, discard_below), sep="", collapse=",")
  print(commands)
  res=eval(parse(text=paste("gcluster.run(",commands,")")))
  res
}


##############
# Create pooled tracks
create_pooled_tracks <- function()
{
  if (sch_haploid) {
    create_pooled_hap_tracks()
  }
  else {
    create_pooled_hyb_tracks()
  }
  
}

####
create_pooled_hap_tracks <- function()
{

  all_nxt_dirs = gsub("_", ".", gsub("scell.nextera.", "", grep("NXT", rownames(sch_chrom_stat), value=T)), fixed=T)
  good_nxt_dirs = gsub("_", ".", gsub("scell.nextera.", "", grep("NXT", sch_good_cells, value=T)), fixed=T)
  all_nst_dirs = gsub("_", ".", gsub("scell.nextera.", "", grep("NST", rownames(sch_chrom_stat), value=T)), fixed=T)
  good_nst_dirs = gsub("_", ".", gsub("scell.nextera.", "", grep("NST", sch_good_cells, value=T)), fixed=T)

  all_dirs = c(paste("2i", all_nxt_dirs, sep="/"), paste("serum", all_nst_dirs, sep="/"))
  good_dirs = c(paste("2i", good_nxt_dirs, sep="/"), paste("serum", good_nst_dirs, sep="/"))

  message("Pooling all...")
  sch_create_pooled_track(all_dirs, pool_track_nm="scell.nextera.pool_all_hap_2i_serum_es", adj_base_dir=sch_base_dir)

  message("Pooling good...")
  sch_create_pooled_track(good_dirs, pool_track_nm="scell.nextera.pool_good_hap_2i_serum_es", adj_base_dir=sch_base_dir)

  message("Pooling all 2i...")
  sch_create_pooled_track(all_nxt_dirs, pool_track_nm="scell.nextera.pool_all_hap_2i_es", adj_base_dir=paste(sch_base_dir, "2i", sep="/"))

  message("Pooling good 2i...")
  sch_create_pooled_track(good_nxt_dirs, pool_track_nm="scell.nextera.pool_good_hap_2i_es", adj_base_dir=paste(sch_base_dir, "2i", sep="/"))

  message("Pooling all serum...")
  sch_create_pooled_track(all_nst_dirs, pool_track_nm="scell.nextera.pool_all_hap_serum_es", adj_base_dir=paste(sch_base_dir, "serum", sep="/"))

  message("Pooling good serum...")
  sch_create_pooled_track(good_nst_dirs, pool_track_nm="scell.nextera.pool_good_hap_serum_es", adj_base_dir=paste(sch_base_dir, "serum", sep="/"))

}

####
create_pooled_hyb_tracks <- function()
{
  create_pooled_hyb_tracks_helper(all_nm="scell.nextera.pool_all_hyb_2i_all_es", good_nm="scell.nextera.pool_good_hyb_2i_all_es", pattern="1CDU|1CDX|1CDES")

  create_pooled_hyb_tracks(all_nm="scell.nextera.pool_all_hyb_2i_no_sort_es", good_nm="scell.nextera.pool_good_hyb_2i_no_sort_es", pattern="1CDES")

  create_pooled_hyb_tracks(all_nm="scell.nextera.pool_all_hyb_2i_sort_es", good_nm="scell.nextera.pool_good_hyb_2i_sort_es", pattern="1CDU|1CDX")

  create_pooled_hyb_tracks(all_nm="scell.nextera.pool_all_hyb_serum_es", good_nm="scell.nextera.pool_good_hyb_serum_es", pattern="1CDS")
}

##
create_pooled_hyb_tracks_helper <- function(all_nm="scell.nextera.pool_all_hyb_es", good_nm="scell.nextera.pool_good_hyb_es", pattern="1CD")
{
  all_names = gsub("_", ".", gsub("scell.nextera.", "", grep("1CD", rownames(sch_chrom_stat), value=T)), fixed=T)
  good_names = gsub("_", ".", gsub("scell.nextera.", "", grep("1CD", sch_good_cells, value=T)), fixed=T)

  all_names = gsub("1CDES.p", "1CDES_p", all_names)
  good_names = gsub("1CDES.p", "1CDES_p", good_names)
  
  all_names = all_names[grep(pattern, all_names, perl=T)]
  good_names = good_names[grep(pattern, good_names, perl=T)]
  
  all_hyb_es_dirs  = sapply(sch_base_dir, paste, all_names,  sep="/")
  good_hyb_es_dirs = sapply(sch_base_dir, paste, good_names, sep="/")

  all_hyb_es_dirs = all_hyb_es_dirs[file.exists(all_hyb_es_dirs)]
  good_hyb_es_dirs = good_hyb_es_dirs[file.exists(good_hyb_es_dirs)]


  message("Pooling all hyb ES...")
  sch_create_pooled_track(all_hyb_es_dirs, pool_track_nm=all_nm, adj_base_dir=NULL)

  message("Pooling good hyb ES...")
  sch_create_pooled_track(good_hyb_es_dirs, pool_track_nm=good_nm, adj_base_dir=NULL)

}

###########
import_rna_seq_fastq <- function(base_dir="/net/mraid14/export/tgdata/db/tgdb/mm9/rawdata/scell_hic/rna_seq/es_hap_and_129Cast", fastq_prefs=c("lane4696_GCCAAT_129_Cast-2_L001_R", "lane4696_ATGTCA_129_Cast-3_L001_R", "lane4696_CGATGT_H129-1_L001_R", "lane4696_GTGAAA_H129-3_L001_R"), track_nms=c("rna.129_Cast.ES.pf_rep1", "rna.129_Cast.ES.pf_rep2", "rna.129_hap.ES.pf_rep1", "rna.129_hap.ES.pf_rep2"), mates_r=c(1,3), binsize=50, min_qual=38, rebuild=F)
{
  odir = paste0(base_dir, "/processed")
  system(paste("mkdir -p", odir))
  
  for (i in seq_along(fastq_prefs)) {
    if (!gtrack.exists(tracks_nms[i]) | rebuild) {
      cmd = sprintf("bowtie2 --local -p 8 -x %s/mm9 -1 %s/%s%d.fastq -2 %s/%s%d.fastq -S %s/%s.sam", Sys.getenv("BOWTIE2_IDX"), base_dir, fastq_prefs[i], mates_r[1], base_dir, fastq_prefs[i], mates_r[2], odir, fastq_prefs[i])
      res = system(cmd)
      system(sprintf("perl -e 'print \"chrom\tstart\tend\tvalue\n\"; ' > %s/%s.tab", odir, fastq_prefs[i]))
      
      system(sprintf("cat %s/%s.sam | perl -ane 'if ( $F[1] & 0x2  and $F[1] & 0x40 and $F[4] >= %d ) { $s = $F[3] < $F[7] ? $F[3] : $F[7]; $e = $s + abs ($F[7]); print \"$F[9]\t$F[2]\t$s\t+\n\"; }' >> %s/%s.mapped_reads", odir, fastq_prefs[i], min_qual, odir, fastq_prefs[i]))
      message("importing ", track_nms[i])
      dirs = strsplit(track_nms[i], ".", fixed=T)[[1]]
      pdir = dirs[1]
      
      gdir.create(pdir, F)
      for (cdir in dirs[c(-1, -length(dirs))]) {
        pdir = paste(pdir, cdir, sep="/")
        gdir.create(pdir, F)
      }
      
      gtrack.import_mappedseq(track_nms[i], "Fraser lab mESC RNA-seq", file=sprintf("%s/%s.mapped_reads", odir, fastq_prefs[i]), binsize=binsize, pileup=150, cols.order=1:4, remove.dups=F)
    }
  }
}

########################
create_ctcf_motif_tracks <- function(rebuild=F, base_dir="motifs", tns=c("CTCF_forward", "CTCF_reverse"), keys=1:2, binsize=1)
{
  
  gdir.create(base_dir, showWarnings=F)

  for (i in seq_along(tns)) {
    otn = sprintf("%s.%s", base_dir, tns[i])
    if (rebuild) {
      gtrack.rm(otn, T)
    }
    if (!gtrack.exists(otn)) {
      message(sprintf("Creating %s track...", otn))
      gtrack.create_pwm_energy(otn, "CTCF pwm track", 'ctcf', keys[i], 0, binsize)
    }
  }

}

#######
#
gen_post_mitotic_insu_and_shuffle <- function(mit_tn=paste0(pool_tn, "_group_1_post_m"), min_diag_d=8192, res=2e4, ignore_below=1024)
{
  new_track=sprintf("%s_ins_discard%d_%ds", mit_tn, min_diag_d, ins_scale)
  if (gtrack.exists(new_track)) {
    gtrack.rm(new_track, T)
  }
  gtrack.2d.gen_insu_track(mit_tn, ins_scale, res=res, min_diag_d=min_diag_d, new_track=new_track, ignore_below=ignore_below)

  options(shaman.sge_support=1)
  shaman_shuffle_hic_track(track_db=sch_groot, obs_track_nm=mit_tn, work_dir=sprintf("%s/", sch_rdata_dir))

  mit_s_tn = sprintf("%s_shuffle", mit_tn)

  new_track=sprintf("%s_ins_discard%d_%ds", mit_s_tn, min_diag_d, ins_scale)
  if (gtrack.exists(new_track)) {
    gtrack.rm(new_track, T)
  }
  gtrack.2d.gen_insu_track(mit_s_tn, ins_scale, res=res, min_diag_d=min_diag_d, new_track=new_track, ignore_below)
  
}



