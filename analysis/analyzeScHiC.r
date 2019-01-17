# sc HiC analysis scripts
#
# this assumes one analysis workspace, with global variables showing up with
# the prefix sch_. We balance between native misha calls to tracks and
# manipulation of cached tables with statistics
#
#
# variables:
#
# sch_glob_decay - storing coverage per 1D dist per cell
# this can be used for global filtering, deciding which cells to use
#
# sch_chrom_stat - storing coarse grained decay per chromosome per cell
# to identify cells with missing chromosome and do general stat on
# correlation between different chromosome stats
#
# sch_chrom_marg, sch_chrom_marg_z - total chrom cov per cell, or z-score
# of biase from the mean
#
# sch_filt_cells - list of cells to be used in downstream analysis
#
#
#
# commands: (TBD: finish documenting)
#
# sch_load_config..
#
#
#
library(misha)
#library(Matrix) 
#library(gtools)
library(plyr)
library(dplyr)
library(tidyr)
library(KernSmooth)
library(RColorBrewer)

source("../map3c/TG3C/analyzeHiC.r")
source("utility_fxns.r")
source("insulation.r")
source("ordering.r")

#===================================

# spec_params_fn needs to be set to path to the params file. e.g.,
# spec_params_fn_both = "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/serum_2i/serum_2i_params.r"

sch_load_config <- function(spec_params_fn)
{
    #track base name
    sch_track_base <<- "scell.nextera.NXT_[0-9]"

    # batch file name
    sch_batch_fn <<- "cell_batch_no_EGO.txt"

    #sch tor track
    sch_tor_track <<- "Encode.esd3.replichip.rep1"

    #table directory
    sch_table_dir <<- "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/tables"
    sch_fig_dir <<- "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/figs"
    sch_rdata_dir <<- "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/rdata"
    # haploid or diploid
    sch_haploid <<- T

    #downsamp param
    sch_min_cov <<- 10000

    #max contacts
    sch_max_cov <<- 1e6

    #max self level
    sch_max_non_digest <<- 0.5

    #max trans contacts
    sch_max_trans <<- 0.2

    #max_chromosome coverage bias (for filtering)
    sch_karyo_z_thresh <<- 60

    #max_chromosome coverage bias (for filtering)
    sch_karyo_log2_enr_thresh <<- 1

    #number of clusters for kmeans by decay
    sch_decay_k <<- 13

    # data dir (misha db, redb, external files)
    sch_data_dir <<- "/net/mraid14/export/data/db/tgdb/mm9/"
    
    # misha db dir
    sch_groot <<- sprintf("%s/trackdb", sch_data_dir)

    # redb dir
    sch_redb_dir <<- sprintf("%s/redb", sch_data_dir)

    # external files dir
    sch_extfiles_dir <<- sprintf("%s/rawdata", sch_data_dir)

    # colors
    twenty_colors <<- c("cadetblue", "lightsalmon", "mediumseagreen", "cornflowerblue","bisque3",
            "lightblue3","hotpink3", "steelblue", "indianred1", "aquamarine4",
            "mediumpurple","darkorange2", "palegreen3", "burlywood4", "lightpink3",
            "skyblue4", "darkgoldenrod2", "paleturquoise3", "tomato","thistle3")

    qualitative_colors <<- c("#8dd3c7", "#1f78b4","#ffd92f", "#fb9a99","#33a02c","#cab2d6","#e31a1c","#a6cee3","#ff7f00","#ffff99","#6a3d9a","#b2df8a","#b15928", "#e7298a", "#fdbf6f", "#1b9e77", "#e78ac3", "blue4", "orangered", "violetred4")
    mini_qualitative_colors <<- c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00','#ffff33')

    nettas_trio <<- c("#87FFFF", "black", "#FF413D")
    nettas_trio2 <<- c("#87FFFF", "white", "#FF413D")

    seq_colspec <<- brewer.pal(n=11, "RdYlBu")[7:11]
    diver_colspec <<- rev(brewer.pal(n=11, "RdYlBu"))
    diver_dark_mid_colspec <<- rev(brewer.pal(n=11, "RdYlBu"))[c(1:4, 8:11)]

    darkblue_col <<- seq_colspec[ length(seq_colspec) ]
    
# domains
    pool_tn <<- "scell.nextera.pool_all"
    ins_scale <<- 3e5
    ins_dom_thresh <<- 0.1
    domain_size_range <<- c(2e4, 4e6)
    intra_domain_na_segment_len_to_filter <<- 25e3

    extend_border_region <<- 1e4
    ins_tn <<- sprintf("%s_ins_%ds", pool_tn, ins_scale)

# sch active chromatin track
    active_tn <<- "Encode.esb4.h3k4me3.rep1"
    lamin_tn <<- "laminB1.ESC"

    # insulation analysis
    insu_ds_quant <<- 0.5
    insu_res <<- 1e5
    insu_scale <<- 2e6
    insu_filter_by <<- c("X2M", "far")

    # min contacts to downsample per cell in domains marginal coverage analysis
    sch_dom_cov_ds <<- 5e4
    sch_dom_cov_ds_q <<- 0.2

    # cis decay - high res analysis
    sch_decay_step <<- 0.125
    sch_cis_breaks <<- c(0, 2**seq(10, 27.75, by=sch_decay_step))

    # contact distance grouping thresholds
    sch_remove_near_dist <<- 2**14.5
    sch_far_dists <<- 2**c(22.125, 27.75)
    sch_near_dists <<- 2**c(14.5, 22)
    sch_domain_dists <<- 2**c(16.25, 21)
    #sch_low_near_dists <<- 2**c(12.5, 15)
    sch_high_near_dists <<- 2**c(15, 22)

    # grouping thresholds
    sch_phasing_criteria <<- "newer"

    sch_glob_min_wm_dist_bin <<- 2 ** 15.5
    sch_high_near_s_start <<- 0.7
    sch_high_near_mid_s   <<- 0.82
    sch_mitotic_cutoff    <<- 0.015
    sch_mitotic_bin       <<- 2**20.875

    # global batches specific params
    sch_cells_group <<- "hap_serum_2i"

    # trans mega covered regions binsize
    sch_megacov_iter <<- 5e4

    # chroms to use
    sch_chrom_levels <<- c(1:19, 'X', 'Y', 'M')
    sch_chroms <<- c(1:19, 'X')
    
    if (!is.null(spec_params_fn)) {
      source(spec_params_fn)
    }

    sch_use_batches <<-  switch(sch_cells_group,
                                hap_serum_2i=c("1", "1268", "1527", "1632", "1654", "1749", "2434", "2435", "1835", "1845", "1994", "2005", "2014", "2027", "3462", "3463", "3580", "3751", "3772", "3941"),
                                hyb_2i_idx_sort_es=c("1CDU_2312", "1CDU_2851", "1CDU_2863", "1CDU_2891", "1CDX1_3136", "1CDX2_3136", "1CDX3_3136", "1CDX4_3136", "1CDX1_3158", "1CDX1_3265", "1CDX2_3197", "1CDX2_3198", "1CDX3_3221", "1CDX3_3222", "1CDX4_3266", "1CDX4_3267", "1CDES_p1", "1CDES_p3", "1CDES_p4", "1CDES_p5", "1CDES_p6", "1CDES_p7", "1CDES_p8", "1CDES_p9", "1CDES_p10"),
                                hyb_2i_idx_sort_serum_es=c("1CDS1_2812", "1CDS1_2891", "1CDS1_2948", "1CDS1_2950", "1CDS1_3046", "1CDS2_2812", "1CDS2_2947", "1CDS2_2949", "1CDS2_2951", "1CDS2_3046", "1CDU_2312", "1CDU_2851", "1CDU_2863", "1CDU_2891", "1CDX1_3136", "1CDX2_3136", "1CDX3_3136", "1CDX4_3136", "1CDX1_3158", "1CDX1_3265", "1CDX2_3197", "1CDX2_3198", "1CDX3_3221", "1CDX3_3222", "1CDX4_3266", "1CDX4_3267", "1CDES_p1", "1CDES_p3", "1CDES_p4", "1CDES_p5", "1CDES_p6", "1CDES_p7", "1CDES_p8", "1CDES_p9", "1CDES_p10"))


    sch_cond_colors <<- switch(sch_cells_group,
                               hap_serum_2i=list("2i_all"="red", "2i_G1"="dodgerblue1", "Serum_G1/S"="orange"),
                               hyb_2i_idx_sort_es=list("1CDU"="red", "1CD_G1"="blue", "1CD_eS"=rgb(0.7, 0.3, 0.7), "1CD_mS"=rgb(0.48, 0.76, 0.48), "1CD_lS_G2"=rgb(0, 0.37, 0), "NoSort"='darkgray'),
                               hyb_2i_idx_sort_serum_es=list("1CDU"="red", "1CDS1"="dodgerblue1", "1CDS2"="orange", "1CD_G1"="blue", "1CD_eS"=rgb(0.7, 0.3, 0.7), "1CD_mS"=rgb(0.48, 0.76, 0.48), "1CD_lS_G2"=rgb(0, 0.37, 0), "NoSort"='darkgray')
                               )

    # misha options
    gdb.init(sch_groot)
    options(gmax.data.size=1e8)

    .update_shaman_options()
}

#===================================

sch_build_all <- function(spec_params_fn)
{
    sch_load_config(spec_params_fn)

    system(sprintf("mkdir -p %s", sch_table_dir))
    system(sprintf("mkdir -p %s", sch_fig_dir))
    system(sprintf("mkdir -p %s", sch_rdata_dir))
    
    sch_create_glob_decay()
    sch_create_glob_decay_res()
    sch_create_repli_stat()
    sch_create_chrom_stat_marg_and_z()
    sch_create_fend_dup_stat()
    count_dup_trans_contacts_across_all_cells()
    compute_chrom_pairs_total_contacts()
    count_unique_fends_per_cell()
    compute_reads_per_contact()
    
    sch_load_tables(spec_params_fn)

    sch_report_karyo_and_qc()
    sch_report_tor_chromcov()
    sch_report_dup_stats(n1_b=NULL)
}

#===================================

sch_load_tables <- function(spec_params_fn)
{
    sch_load_config(spec_params_fn)

    sch_glob_decay <<- read.table(sprintf("%s/glob_decay.txt", sch_table_dir),header=T)
    sch_glob_decay_res <<- read.table(sprintf("%s/glob_decay_res_step%s.txt", sch_table_dir, as.character(sch_decay_step)),header=T, sep="\t")
    sch_chrom_decay_res <<- read.table(sprintf("%s/chrom_decay_res_step%s_p1_of1.txt", sch_table_dir, as.character(sch_decay_step)),header=T, sep="\t")
    sch_chrom_stat <<- read.table(sprintf("%s/chrom_stat.txt", sch_table_dir),header=T)
    sch_chrom_marg <<- read.table(sprintf("%s/chrom_marg.txt", sch_table_dir),header=T)
    sch_chrom_marg_z <<- read.table(sprintf("%s/chrom_marg_z.txt", sch_table_dir),header=T)
    sch_chrom_marg_enr <<- read.table(sprintf("%s/chrom_marg_enr.txt", sch_table_dir),header=T)
    sch_tor_stat <<- read.table(sprintf("%s/tor_stat.txt", sch_table_dir),header=T)
    sch_fend_dup <<- read.table(sprintf("%s/fend_cis_dup.txt", sch_table_dir), header=T)

    sch_trans_switchers <<- load_table_if_exists(sprintf("%s/trans_dup_contacts.txt", sch_table_dir), header=T)
    sch_trans_pairs_contact_count <<- load_table_if_exists(sprintf("%s/sch_trans_pairs_contact_count_mul.txt", sch_table_dir), header=T)
    sch_domain_raw_marg_cov_ds <<- load_table_if_exists(sprintf("%s/domain_raw_marg_cov_ds%.2f.txt", sch_table_dir, 0.1), header=T)
    sch_domain_raw_marg_cov <<-  load_table_if_exists(sprintf("%s/domain_raw_marg_cov.txt", sch_table_dir), header=T)

    sch_cis_contact_mul <<- load_table_if_exists(sprintf("%s/sch_cis_contact_mul.txt", sch_table_dir), header=T)
    sch_trans_contact_mul <<- load_table_if_exists(sprintf("%s/sch_trans_contact_mul.txt", sch_table_dir), header=T)
    sch_batch <<- read.table(sprintf("%s/%s", sch_table_dir, sch_batch_fn),header=T, stringsAsFactors=F)
    ub = unique(sch_batch)
    rownames(ub) = ub$batch
    sch_batch_colors <<- unlist(sch_cond_colors[ub[1:max(ub$batch), 'cond']])

    # calculate (and write to file) tor stats on the fly because they rely on
    # sch_good_cells - want to preserve flexibility to modify qc thresholds on the fly
    sch_select_qc_cells()
    sch_compute_tor_marg()
    
    # filter by decay analysis
    sch_decay_metrics <<- group_cells_by_metrics(plot_extra=F, criteria=sch_phasing_criteria, cells=sch_good_cells)

    sch_good_cells_ordered <<- intersect(rownames(sch_decay_metrics), sch_good_cells)
    sch_good_cells <<- intersect(sch_good_cells, rownames(sch_decay_metrics))
    
    # load ab compartment if exist
    sch_ab_counts <<- load_table_if_exists(sprintf("%s/cell_ab_counts.txt", sch_table_dir))

}

#===================================

.sch_gen_coord_chrom_base <- function(intervs, bsz)
{
    cmax = tapply(intervs$end, as.character(intervs$chrom), max)
    cmin = tapply(intervs$start, as.character(intervs$chrom), min)
    chrnum = as.numeric(sub("M", "22", sub("Y", "21", sub("X", "20", sub("chr", "", names(cmax))))))

    cdelta = cmax[order(chrnum)]-cmin[order(chrnum)]
    cdelta = ceiling(cdelta/bsz)
    base = 1+cumsum(c(0,cdelta[-length(cdelta)]))
    names(base) = names(cdelta)

    return(list(chr_base = base, max_id = sum(cdelta)))
}

#===================================

sch_create_glob_decay <- function()
{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   CREATES     sch_glob_decay
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    message("sch_create_glob_decay")
    nms = gtrack.ls(sch_track_base)
    stats = data.frame(c0=c(), c1 = c(), c2 =c(), c3 =c(), c4 =c(), c5 = c(), ctrans = c())

    commands = paste0(".sch_gen_mat_decay(nm=\"", nms, "\")", collapse=",")

    res = eval(parse(text=paste("gcluster.run(",commands,")")))

    for (i in seq_along(res)) {
      rv = res[[i]]$retv
      stopifnot(!is.null(rv), typeof(rv) == "list")
      stats = rbind(stats, rv)
    }

    rownames(stats) = nms

    sch_glob_decay <<- stats

    write.table(sch_glob_decay, sprintf("%s/glob_decay.txt", sch_table_dir),quote=F, sep="\t")
}


#===================================

sch_create_glob_decay_res <- function(n_partitions=1, chroms=sch_chroms)
{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   CREATES     sch_glob_decay_res
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    message("sch_create_glob_decay_res")
    nms = gtrack.ls(sch_track_base)


    commands = paste0(".sch_gen_mat_decay_res(nm=\"", nms, "\", chroms=chroms, n_partitions=n_partitions)", collapse=",")

    res = eval(parse(text=paste("gcluster.run(",commands,")")))
    vals = c()
    for (i in seq_along(res)) {
      rv = res[[i]]$retv
      stopifnot(!is.null(rv), typeof(rv) == "integer" | typeof(rv) == "double")
      vals = c(vals, as.vector(rv))
    }

    stats = array(vals,  dim=c(n_partitions, length(chroms) * length(sch_cis_breaks), length(nms)), dimnames=list(part=1:n_partitions, chrom_bin=apply(expand.grid(paste0("chr", chroms), c(sch_cis_breaks[-1], "trans")), 1, function(v) paste(v, collapse="_")), nm=nms))
    

    sch_chrom_decay_res <<- apply(stats, c(3,2), sum)
    
    y <<- t(apply(sch_chrom_decay_res, 1, tapply, gsub("chr.*_", "X", colnames(sch_chrom_decay_res), perl=T), sum))
    colnames(y)[grep("trans", colnames(y))] = "trans"
    
    #rownames(y) = nms
    cols = c(paste0('X', as.character(sch_cis_breaks[-1])), "trans")
    sch_glob_decay_res <<- y[, cols]

    # when summing genome-wide - fix double counting of trans contacts
    sch_glob_decay_res[, 'trans'] = sch_glob_decay_res[, 'trans'] / 2
    
    for (i in 1:n_partitions) {
      write.table(t(stats[i, , ]), sprintf("%s/chrom_decay_res_step%s_p%d_of%d.txt", sch_table_dir, as.character(sch_decay_step), i, n_partitions), quote=F, sep="\t")      
    }

    write.table(sch_glob_decay_res, sprintf("%s/glob_decay_res_step%s.txt", sch_table_dir, as.character(sch_decay_step)), quote=F, sep="\t")


}


#===================================

sch_create_chrom_stat_marg_and_z <- function()
{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   CREATES     sch_chrom_stat
    #               sch_chrom_marg
    #               sch_chrom_marg_z
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    message("sch_create_chrom_stat_marg_and_z")
    nms = gtrack.ls(sch_track_base)

    #chr_counts = data.frame()
    #for (nm in nms) {
    #    message("collect ", nm)
    #    ccount = .sch_gen_mat_chr_cov(nm)
    #    chr_counts = rbind(chr_counts, ccount)
    #}

    commands = paste0(".sch_gen_mat_chr_cov(nm=\"", nms, "\")", collapse=",")
    res = eval(parse(text=paste("gcluster.run(",commands,")")))
    chr_counts = do.call("rbind", lapply(res, function(x) { x$retv } ))

    n_chroms = nrow(gintervals.all())
    rownames(chr_counts) = nms
    colnames(chr_counts) = paste(c(rep("self",n_chroms),rep("2M", n_chroms),
                rep("far", n_chroms),rep("trans", n_chroms)),
                rep(levels(gintervals.all()$chrom),4),
                sep="_")

    # calculate "sch_chrom_stat" -- number of contacts per chromosome in each
    # band named above
    
    chr_counts = chr_counts[, grep("chrM|chrY", colnames(chr_counts), perl=T, invert=T)]
    
    sch_chrom_stat <<- chr_counts
    write.table(sch_chrom_stat, sprintf("%s/chrom_stat.txt", sch_table_dir),quote=F, sep="\t")

    # calculate "sch_chrom_marg" -- total number of contacts per chromosome, and
    # "sch_chrom_marg_z" -- a z score of sorts for the total coverage (comparing
    # observed to population expected counts per chromosome)
    marg_and_z <- sch_calc_chrom_marg_z(band=c("total"))
    sch_chrom_marg <<- marg_and_z[["marg"]]
    sch_chrom_marg_z <<- marg_and_z[["marg_z"]]
    sch_chrom_marg_enr <<- marg_and_z[["enr"]]
    write.table(sch_chrom_marg, sprintf("%s/chrom_marg.txt", sch_table_dir),quote=F, sep="\t")
    write.table(sch_chrom_marg_z, sprintf("%s/chrom_marg_z.txt", sch_table_dir),quote=F, sep="\t")
    write.table(sch_chrom_marg_enr, sprintf("%s/chrom_marg_enr.txt", sch_table_dir),quote=F, sep="\t")

}

#===================================

sch_calc_chrom_marg_z <- function(band=c("total"), cells=gtrack.ls(sch_track_base))
{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   REQUIRES    sch_chrom_stat
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    # the arg <band> is either a string or a vector of strings
    # vector can contain the vals "self", "2M", "far", "trans" or "total"

    if (band[1]=="total") {
        band <- c("self", "2M", "far", "trans")
    }

    # filter chrom stat matrix to only the cells of interest
    filt_chrom_stat <- sch_chrom_stat[cells, ]

    # prealloc chrom_marg matrix
    chrom_marg <- matrix(0, nrow=dim(filt_chrom_stat)[1], ncol=dim(filt_chrom_stat)[2]/4)

    for (this_band in band) {
        chrom_marg <- chrom_marg + filt_chrom_stat[,grep(this_band, colnames(filt_chrom_stat))]
    }


    colnames(chrom_marg) <- sub("^X", "", sub(paste0(band[1], "_"), "", colnames(chrom_marg)))
    chrom_p = colSums(chrom_marg)/sum(chrom_marg)
    chrom_exp = rowSums(chrom_marg) %*% t(chrom_p)
    chrom_marg_z <- (chrom_marg - chrom_exp) / sqrt(chrom_marg)
    chrom_enr <- log2( (chrom_marg + 1) / (chrom_exp + 1) )

    return(list(marg=as.matrix(chrom_marg), marg_z=as.matrix(chrom_marg_z), enr=as.matrix(chrom_enr)))
}

#===================================

sch_create_fend_dup_stat <- function()
{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   CREATES     sch_fend_dup
    #
    #   Fraction of fends with > 1 cis contacts out of covered fends per chromosome
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    message("sch_create_fend_dup_stat")
    nms = gtrack.ls(sch_track_base)

    ints = gintervals.2d.all()
    ints = ints[as.character(ints$chrom1) == as.character(ints$chrom2), ]
    ints = ints[grep("chrM|chrY", ints$chrom1, invert=T), ]

    stats = matrix(0, length(nms), nrow(ints))
    rownames(stats) = nms
    colnames(stats) = as.character(ints$chrom1)

    commands = paste0(".sch_create_fend_dup_nm(nm=\"", nms, "\", ints=ints)", collapse=",")

    res = eval(parse(text=paste("gcluster.run(",commands,")")))

    for (i in seq_along(res)) {
      rv = res[[i]]$retv
      stopifnot(!is.null(rv), typeof(rv) == "double")
      stats[nms[i], names(rv)] = rv
    }

    sch_fend_dup <<- stats

    write.table(sch_fend_dup, sprintf("%s/fend_cis_dup.txt", sch_table_dir),quote=F, sep="\t")
}

.sch_create_fend_dup_nm <- function(nm, ints)
{
  x = gextract(nm, intervals=ints)
  tapply(x$start1, as.character(x$chrom1), function(y) { z = table(y); sum(z > 1) / length(z) })
}

#===================================

.sch_gen_mat_decay <- function(nm)
{
    mat = gextract(nm, gintervals.2d.all())
    #cis rate
    f_cis = mat$chrom1==mat$chrom2

    #10M-40/2.5-10
    mat$d = ifelse(f_cis, abs(mat$start1-mat$start2), 1e+9)

    ctrans = sum(mat$d == 1e+9)
    c5 = sum(mat$d>1e+7 & mat$d<1e+8)
    c4 = sum(mat$d>1e+6 & mat$d<1e+7)
    c3 = sum(mat$d>1e+5 & mat$d<1e+6)
    c2 = sum(mat$d>1e+4 & mat$d<1e+5)
    c1 = sum(mat$d>1e+3 & mat$d<1e+4)
    c0 = sum(mat$d>0 & mat$d<1e+3)

    return(list(c0=c0, c1=c1,c2=c2,c3=c3,c4=c4,c5=c5,ctrans=ctrans))
}

#===================================

.sch_gen_mat_decay_res <- function(nm, chroms=sch_chroms, n_partitions=1)
{
  col_names = apply(expand.grid(paste0("chr", chroms), c(sch_cis_breaks[-1], "trans")), 1, function(v) paste(v, collapse="_"))

  
  va = matrix(0, n_partitions, length(col_names))
  colnames(va) = col_names

  x = partition_contacts(nm, n_partitions, chroms)
  x = x[ x$chrom1 != x$chrom2 | x$start2 > x$start1, ]

  if (!is.null(x)) {
    x$dist_bin = cut(ifelse(x$chrom1 == x$chrom2, x$start2 - x$start1, 1e9), breaks=c(sch_cis_breaks, 1e9), labels=c(sch_cis_breaks[-1], "trans"))
    x$bin = paste(x$chrom1, x$dist_bin, sep="_")
    
    v = table(x$part, x$bin)

    va[as.numeric(rownames(v)), colnames(v)] = va[as.numeric(rownames(v)), colnames(v)]  + v
  }

  return(va)
}

#===================================

.sch_gen_mat_chr_cov <- function(nm)
{
    mat = gextract(nm, gintervals.2d.all())
    #cis rate
    f_cis = mat$chrom1==mat$chrom2

    #10M-40/2.5-10
    mat$d = ifelse(f_cis, abs(mat$start1-mat$start2), 1e+9)
    mat$mode = ifelse(mat$chrom1==mat$chrom2, ifelse(mat$d<1000, 0, ifelse(mat$d<2e+6, 1, 2)), 3)

    dummy_counts = expand.grid(gintervals.all()$chrom, 0:3)
    chr_counts = as.vector(table(c(mat$chrom1, dummy_counts[,1]), c(mat$mode, dummy_counts[,2])))

    return(chr_counts - 1)
}

#===================================

sch_create_repli_stat <- function(n_partitions=10, chroms=sch_chroms)
{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   CREATES     sch_tor_stat
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  message("sch_create_repli_stat")
  n_chroms = length(chroms)
  nms = gtrack.ls(sch_track_base)

  tor_counts = array(0, dim=c(n_partitions, length(nms), n_chroms * 8), dimnames=list(part=1:10, nm=nms, mode=paste(c(rep("self_tor1",n_chroms),rep("self_tor2", n_chroms),
                rep("self_tor3",n_chroms),rep("self_tor4", n_chroms),
                rep("good_tor1", n_chroms),rep("good_tor2", n_chroms),
                rep("good_tor3", n_chroms),rep("good_tor4", n_chroms)),
                rep(paste0("chr", chroms),8),
                sep="_")))
  
  commands = paste0(".sch_create_repli_nm(nm=\"", nms, "\", chroms=chroms, n_partitions=n_partitions)", collapse=",")

  res = eval(parse(text=paste("gcluster.run(",commands,")")))

  for (i in seq_along(res)) {
    cov = res[[i]]$retv
    if (!is.null(cov)) {
      stopifnot(typeof(cov) == "list")
      
      tor_counts[ cbind(cov$Var1, nms[i], paste(cov$Var3, cov$Var2, sep="_"))] = cov$Freq
    }
  }

  
  sch_tor_stat <<- apply(tor_counts, 2:3, sum)

  write.table(sch_tor_stat, sprintf("%s/tor_stat.txt", sch_table_dir),quote=F, sep="\t")
  
  for (i in 1:n_partitions) {
    write.table(tor_counts[i, , ], sprintf("%s/tor_stat_p%d_of%d.txt", sch_table_dir, i, n_partitions), quote=F, sep="\t")      
  }

}

.sch_create_repli_nm <- function(nm, n_partitions=1, chroms=sch_chroms)
{
  gvtrack.create("tor_proj", sch_tor_track)
  gvtrack.iterator("tor_proj", sshift=-10000, eshift=10000, dim=1)

  mat = partition_contacts(nm, n_partitions, chroms)

  str_chrs = paste0("chr", chroms)

  scope = gintervals.2d.all()
  scope = scope[is.element(scope$chrom1, str_chrs) & is.element(scope$chrom2, str_chrs), ]
  mat$tor_proj = gextract("tor_proj", intervals=scope, iterator=nm)$tor_proj
  mat = mat[!is.na(mat$tor_proj), ]
  
  cis = mat$chrom1 == mat$chrom2

  mat$mode = floor(pmin(pmax(mat$tor_proj, -2), 2 - 1e-4))
  
  mat$mode = paste0(ifelse(cis & abs(mat$start1-mat$start2) < 1e+3, "self", "good"), "_tor", mat$mode + 3)

  as.data.frame(table(mat$part, as.character(mat$chrom1), mat$mode))
}


#===================================

sch_compute_tor_marg <- function()
{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   REQUIRES    sch_good_cells
    #               sch_tor_stat
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   CREATES     sch_tor_marg
    #               sch_tor_marg_n
    #               sch_nodig_tor_marg
    #               sch_nodig_tor_marg_n
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    tor1_cols = grep("tor1", names(sch_tor_stat))
    tor2_cols = grep("tor2", names(sch_tor_stat))
    tor3_cols = grep("tor3", names(sch_tor_stat))
    tor4_cols = grep("tor4", names(sch_tor_stat))

    sch_tor_marg_all <<- data.frame(
            tor1 = rowSums(sch_tor_stat[,tor1_cols]),
            tor2 = rowSums(sch_tor_stat[,tor2_cols]),
            tor3 = rowSums(sch_tor_stat[,tor3_cols]),
            tor4 = rowSums(sch_tor_stat[,tor4_cols]))
    rownames(sch_tor_marg_all) <<- rownames(sch_tor_stat)

    sch_tor_marg_all_n <<- sch_tor_marg_all/rowSums(sch_tor_marg_all)

    sch_tor_marg <<- sch_tor_marg_all[sch_good_cells,]
    sch_tor_marg_n <<- sch_tor_marg/rowSums(sch_tor_marg)

    tor1_self_cols = grep("self_tor1", names(sch_tor_stat))
    tor2_self_cols = grep("self_tor2", names(sch_tor_stat))
    tor3_self_cols = grep("self_tor3", names(sch_tor_stat))
    tor4_self_cols = grep("self_tor4", names(sch_tor_stat))

    sch_nodig_tor_marg <<- data.frame(
            tor1 = rowSums(sch_tor_stat[sch_good_cells,tor1_self_cols]),
            tor2 = rowSums(sch_tor_stat[sch_good_cells,tor2_self_cols]),
            tor3 = rowSums(sch_tor_stat[sch_good_cells,tor3_self_cols]),
            tor4 = rowSums(sch_tor_stat[sch_good_cells,tor4_self_cols]))
    rownames(sch_nodig_tor_marg) <<- sch_good_cells

    sch_nodig_tor_marg_n <<- sch_nodig_tor_marg/rowSums(sch_nodig_tor_marg)

}

#===================================

sch_select_qc_cells <- function()
{
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   REQUIRES    sch_chrom_marg_z
    #               sch_glob_decay
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #   CREATES     sch_good_cells
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #bad_karyo_cells = rownames(sch_chrom_marg_z)[apply(abs(sch_chrom_marg_z), 1, max) > sch_karyo_z_thresh]
    bad_karyo_cells = rownames(sch_chrom_marg_enr)[apply(abs(sch_chrom_marg_enr), 1, max) > sch_karyo_log2_enr_thresh]

    sch_good_cells <<- rownames(sch_glob_decay)[sch_glob_decay$c0/rowSums(sch_glob_decay) < sch_max_non_digest]

    sch_good_cells <<- setdiff(sch_good_cells, bad_karyo_cells)
    sch_good_cells <<- setdiff(sch_good_cells, grep("b$", sch_good_cells, val=T))
    sch_good_cells <<- setdiff(sch_good_cells, rownames(sch_glob_decay)[rowSums(sch_glob_decay)/2 < sch_min_cov])
    sch_good_cells <<- setdiff(sch_good_cells, rownames(sch_glob_decay)[rowSums(sch_glob_decay)/2 > sch_max_cov])
    sch_good_cells <<- setdiff(sch_good_cells, rownames(sch_glob_decay)[sch_glob_decay$ctrans/rowSums(sch_glob_decay)> sch_max_trans])

    if (is.element("valid", colnames(sch_batch))) {
      sch_good_cells <<- setdiff(sch_good_cells, rownames(sch_batch[sch_batch$valid == 0, ]))
    }
    
    #cell_num = as.numeric(sub(sch_track_base,"", sch_good_cells))
    sch_good_cells <<- sch_good_cells[mixedorder(sch_good_cells)]
}


#===================================

sch_report_karyo_and_qc <- function(trans_trim_to=NULL)
{
  cells = rownames(sch_batch)
  cells = intersect(cells, rownames(sch_glob_decay))

  batch_numbers = unique(sch_batch$batch_id)

  ind = c(1, 1+cumsum(rle(sch_batch$batch)$lengths))
  ind = ind[-length(ind)]

  far_close_ratio = log2(sch_glob_decay[,6]/(sch_glob_decay[2] + sch_glob_decay[3]))
  png(sprintf("%s/qc_far_close_ratio.png", sch_fig_dir), width=1500, height=600)
  par(mar=c(7,6,6,0))
  plot(far_close_ratio[cells,],
       pch=21, bg=sch_batch_colors[sch_batch[cells,1]], col="black",
       xaxt="n",
       xlab="",
       ylab="",
       cex.axis=2)
  abline(h=0, col='black')
                                        #title(xlab="Batch number", cex.lab=2, line=5)
  title(ylab="Far cis (1e7 < d < 1e8 ) / close cis (1e4 < d < 1e6) (log2)", cex.lab=2, line=4)
  abline(v=ind, col="gray")
  axis(side=1, at=ind[c(TRUE, FALSE)], labels=batch_numbers[c(TRUE, FALSE)], las=2, cex.axis=1)
  axis(side=3, at=ind[c(FALSE, TRUE)], labels=batch_numbers[c(FALSE, TRUE)], las=2, cex.axis=1)
  dev.off()

  png(sprintf("%s/qc_tot_cov.png", sch_fig_dir), width=1500, height=600)
                                        # batch_colors = colorRampPalette(c("black", "brown","yellow", "pink", "blue", "purple", "green", "red"))(max(sch_batch$batch))
  par(mar=c(7,6,6,0))
  plot(pmin(rowSums(sch_glob_decay)[cells]/2, sch_max_cov),
       pch=21, bg=sch_batch_colors[sch_batch[cells,1]], col="black",
       xaxt="n",
       xlab="",
       ylab="",
       cex.axis=2)
                                        #title(xlab="Batch number", cex.lab=2, line=5)
  title(ylab="Total coverage", cex.lab=2, line=4)
  abline(v=ind, col="gray")
  axis(side=1, at=ind[c(TRUE, FALSE)], labels=batch_numbers[c(TRUE, FALSE)], las=2, cex.axis=1)
  axis(side=3, at=ind[c(FALSE, TRUE)], labels=batch_numbers[c(FALSE, TRUE)], las=2, cex.axis=1)
  dev.off()

  if (is.null(trans_trim_to)) {
    trans_trim_to = max(sch_glob_decay[cells, 'ctrans']/rowSums(sch_glob_decay[cells,]))
  }

  png(sprintf("%s/qc_trans_frac.png", sch_fig_dir), width=1500, height=600)
  par(mar=c(7,6,6,0))
  plot(pmin(sch_glob_decay[cells, 'ctrans']/rowSums(sch_glob_decay[cells,]), trans_trim_to),
       pch=21, bg=sch_batch_colors[sch_batch[cells,1]], col="black",
       xaxt="n",
       xlab="",
       ylab="",
       cex.axis=2)
                                        #title(xlab="Batch number", cex.lab=2, line=5)
  title(ylab="%trans", cex.lab=2, line=4)
  abline(v=ind, col="gray")
  abline(h=sch_max_trans, lty=2, lwd=2)
  axis(side=1, at=ind[c(TRUE, FALSE)], labels=batch_numbers[c(TRUE, FALSE)], las=2, cex.axis=1)
  axis(side=3, at=ind[c(FALSE, TRUE)], labels=batch_numbers[c(FALSE, TRUE)], las=2, cex.axis=1)
  dev.off()

  png(sprintf("%s/qc_karyo_sort.png", sch_fig_dir), width=1000, height=400)
  par(mar=c(7,6,6,0))
  bad_cells = rownames(sch_chrom_marg_z[order(rowSums(abs(sch_chrom_marg_z))),])
  bad_cells_thresh = rownames(sch_chrom_marg_z)[apply(abs(sch_chrom_marg_z), 1, max)> sch_karyo_z_thresh]

  image(pmin(pmax(as.matrix(sch_chrom_marg_z)[bad_cells,],-200),200), col=colorRampPalette(c("red", "white", "blue"))(100), zlim=c(-200,200), yaxt='n')
  axis(4, at=seq(0, 1, length=ncol(sch_chrom_marg_z)), labels=gsub("chr", "", colnames(sch_chrom_marg_z)), las=2)
  abline(v= 1-length(bad_cells_thresh)/nrow(sch_chrom_marg_z),
         col="black", lwd=2)
  dev.off()

  png(sprintf("%s/qc_karyo_nosort.png", sch_fig_dir), width=1000, height=400)
  bad_cells = rownames(sch_chrom_marg_z[order(rowSums(abs(sch_chrom_marg_z))),])
  image(pmin(pmax(as.matrix(sch_chrom_marg_z)[cells,],-200),200), col=colorRampPalette(c("red", "white", "blue"))(100), zlim=c(-200,200), main=sprintf("filt %d cells", length(bad_cells_thresh)), yaxt='n')
  axis(4, at=seq(0, 1, length=ncol(sch_chrom_marg_z)), labels=gsub("chr", "", colnames(sch_chrom_marg_z)), las=2)
  dev.off()


  png(sprintf("%s/qc_nodigest.png", sch_fig_dir), width=1500, height=600)
  par(mar=c(7,6,6,0))
                                        # batch_colors = colorRampPalette(c("black", "brown","yellow", "pink", "blue", "purple", "green", "red"))(max(sch_batch$batch))
  nodig = sch_glob_decay$c0/rowSums(sch_glob_decay)
  plot(nodig[cells],
       pch=21, bg=sch_batch_colors[sch_batch[cells,1]], col="black",
       xaxt="n",
       xlab="",
       ylab="",
       cex.axis=2)
                                        #title(xlab="Batch number", cex.lab=2, line=5)
  title(ylab="% non-digested", cex.lab=2, line=4)
  abline(h=sch_max_non_digest, lty=2)
  abline(v=ind, col="gray")
  axis(side=1, at=ind[c(TRUE, FALSE)], labels=batch_numbers[c(TRUE, FALSE)], las=2, cex.axis=1)
  axis(side=3, at=ind[c(FALSE, TRUE)], labels=batch_numbers[c(FALSE, TRUE)], las=2, cex.axis=1)
  dev.off()

                                        # good_cells = intersect(good_cells,
                                        #       rownames(sch_chrom_marg_z)[apply(abs(sch_chrom_marg_z),
  sch_select_qc_cells()
  good_cells = sch_good_cells

  message("good cells ", length(good_cells))
  png(sprintf("%s/qc_chrom_nodig.png", sch_fig_dir), width=1200, height=900)
  par(mar=c(7,6,6,0))
  layout(matrix(seq(1,20,1),nrow=4,ncol=5))
  par(mar=c(1,1,1,1))
  self_cols = grep("self", names(sch_chrom_stat))
  m_cols = grep("X2M", names(sch_chrom_stat))
  far_cols = grep("far", names(sch_chrom_stat))
  trans_cols = grep("trans", names(sch_chrom_stat))
  chr_self = sch_chrom_stat[good_cells,self_cols]/
    (sch_chrom_stat[good_cells,self_cols] +
     sch_chrom_stat[good_cells,m_cols] +
     sch_chrom_stat[good_cells,far_cols] +
     sch_chrom_stat[good_cells,trans_cols])

  for(col in self_cols[-1]) {
    plot(chr_self[,self_cols[1]], chr_self[, col], pch=21,
         bg=sch_batch_colors[sch_batch[good_cells,1]],
         col="black",
         xlim=c(0,0.5), ylim=c(0,0.5), cex=0.5)
    grid()
  }
  dev.off()

  png(sprintf("%s/qc_nodig_cell.png", sch_fig_dir), width=1000, height=600)
  dig_ord = order(rowSums(chr_self))
  image(pmin(pmax(as.matrix(chr_self[dig_ord,]),0.1),0.4), col=colorRampPalette(c("red", "white", "blue"))(100), zlim=c(0.1,0.4))
  dev.off()

  png(sprintf("%s/qc_nodig_normclust.png", sch_fig_dir), width=1000, height=600)
  chr_self_n = log2(0.01+chr_self/rowMeans(chr_self))
  chr_self_n2 = t(t(chr_self_n)/colMeans(chr_self_n))

  hc = hclust(dist(chr_self_n), "ward.D2")
  hc_chr = hclust(dist(t(chr_self_n2)), "ward.D2")
  image(pmin(pmax(as.matrix(chr_self_n[hc$order,hc_chr$order]),-1),1), col=colorRampPalette(c("red", "white", "blue"))(100), zlim=c(-1,1))
  dev.off()

                                        # cor of no dig and trans to all other dist bins
  dn = sch_glob_decay_res[sch_good_cells, ]

  dn_n = dn / rowSums(dn[,-1])
  dn_n[,1] = dn
  cdn = cor(dn_n)

  png(sprintf("%s/qc_nodig_trans_decay_cor.png", sch_fig_dir), width=800, height=500)
  par(mar=c(7,4,2,1))
  plot(cdn[nrow(cdn),], pch=19, col='blue', lwd=2, ylim=c(-0.5, 0.5), main="cor to dist bins", xlab="", ylab="pearson", xaxt='n')
  axis(1, c(1, seq(4, nrow(cdn)-1, by=3),nrow(cdn)), c("No dig", sapply(as.numeric(gsub("X", "", colnames(dn)[seq(4, ncol(dn)-1, by=3)])), n2str, 2), "Trans"), las=2, cex.axis=0.8)
  points(cdn[1,], pch=19, col='red')
  grid(col='black')
  legend("bottomright", legend=c("Trans", "No dig"), pch=19, col=c("blue", "red"))
  dev.off()

  png(sprintf("%s/qc_karyo_enr_per_chrom.png", sch_fig_dir), width=480, height=280)
  par(mfrow=c(4,5))
  par(mar=c(2,1,1,1))

  for (i in colnames(sch_chrom_marg_enr)) {
    hist(sch_chrom_marg_enr[,i], 500, xlim=c(-2,2), main=i)
    abline(v=c(-1, 1), col='red')
  }
  dev.off()


}

#===================================

sch_report_tor_chromcov <- function()
{
  sch_select_qc_cells()

  good_cells_batch = sch_batch[sch_good_cells,]
  good_cells_batch = good_cells_batch[order(good_cells_batch$batch), ]

  batch_numbers = unique(good_cells_batch$batch_id)
  ind = c(1, 1+cumsum(rle(good_cells_batch$batch)$lengths))
  ind = ind[-length(ind)]

  #coverage on early replicating for different batches
  png(sprintf("%s/tor_all.png", sch_fig_dir), width=1500, height=600)
  par(mar=c(7,6,6,0))

  early_tot_ratio = rowSums(sch_tor_marg_n[rownames(good_cells_batch), c('tor3', 'tor4')])
  nodig_early_tot_ratio = rowSums(sch_nodig_tor_marg_n[rownames(good_cells_batch), c('tor3', 'tor4')])
  plot(early_tot_ratio, pch=21,
       bg=sch_batch_colors[good_cells_batch$batch],
       col="black",
       xaxt="n",
       xlab="",
       ylab="",
       cex.axis=2)
  #title(xlab="Batch number", cex.lab=2, line=5)
  title(ylab="Early / total ratio", cex.lab=2, line=4)
  abline(v=ind, col="gray")
  abline(h=0.55, col="black")
  axis(side=1, at=ind[c(TRUE, FALSE)], labels=batch_numbers[c(TRUE, FALSE)], las=2, cex.axis=1)
  axis(side=3, at=ind[c(FALSE, TRUE)], labels=batch_numbers[c(FALSE, TRUE)], las=2, cex.axis=1)
  points(nodig_early_tot_ratio, pch=21, bg="gray", col="black", cex=0.4)
  dev.off()

  png(sprintf("%s/tor_good_vs_nodig.png", sch_fig_dir), width=600, height=600)
  par(mar=c(10,10,10,10))
  plot(sch_nodig_tor_marg_n$tor4+sch_nodig_tor_marg_n$tor3,
       sch_tor_marg_n$tor4+sch_tor_marg_n$tor3,
       xlim=c(0.45,0.7),ylim=c(0.45,0.7),
       xlab="TOR non-digested",
       ylab="TOR good",
       cex.lab=2, cex.axis=2)
  grid()
  dev.off()

  #we can make a bigger effort in clustering cells to early/late S-phase. However currently
  #it is not clear how useful this is going to be
  # hc = hclust(dist(as.matrix(sch_tor_marg_n2)), "ward.D2")
  # tor_ord = order(sch_tor_marg_n2[,4]+0.2*sch_tor_marg_n2[,3])
  sch_tor_marg_n2 <<- t(t(sch_tor_marg_n)/colMeans(sch_tor_marg_n))
  # tor_toclust = data.frame(el = log2(sch_tor_marg_n2[,3]+sch_tor_marg_n2[,4])-log2(sch_tor_marg_n2[,1]+sch_tor_marg_n2[,2]),
  #             vearly=5*log2(sch_tor_marg_n2[,4]/sch_tor_marg_n2[,4]),
  #             vlate=log2(sch_tor_marg_n2[,1]/sch_tor_marg_n2[,2]))
  # km = kmeans(tor_toclust, 5, iter.max=100, nstart=100)
  sch_tor_marg_n2 = log2(sch_tor_marg_n2)

  tor_ord = order(sch_tor_marg_n2[,4]+0.4*sch_tor_marg_n2[,3])
  png(sprintf("%s/tor_clust.png", sch_fig_dir), width=1000, height=500)
  image(pmax(pmin(as.matrix(sch_tor_marg_n2[tor_ord,]),0.7),-0.7), col=colorRampPalette(c("blue", "white", "red"))(100), zlim=c(-0.7,0.7))
  dev.off()

  #test correlation of chromosome coverage - on late

  #test correlation of chromosome early/late ratio

}
#===================================
sch_report_dup_stats <- function(n1_b= c(1845, 1994, 2005, 2014, 2027))
{
  dup = sch_fend_dup[grep("EGO", sch_good_cells, invert=T, value=T),]

  # dup vs coverage
  total = rowSums(sch_chrom_stat)/2
  mdup = rowMeans(sch_fend_dup)

  ind = total < sch_max_cov
  total = total[ind]
  mdup = mdup[ind]

  png(sprintf("%s/dup_vs_cov.png", sch_fig_dir), width=300, height=300)
  plot(mdup, total, type='p', pch=19, cex=0.5, main="", xlab="mean Pr[dup]", ylab="#contacts")
  points(mdup[sch_good_cells], total[sch_good_cells], pch=19, cex=0.5, col='red')
  legend("topleft", legend=c("all", "good"), col=c("black", "red"), pch=19)
  max_dup = 0.015
  max_cov = 110000
  abline(h=max_cov)
  abline(v=max_dup)
#  text(x=0, y=0.9*max_cov, pos=4, labels=sprintf("%.0f%%", 100* sum(mdup[sch_good_cells] < max_dup & total[sch_good_cells] < max_cov) / length(mdup[sch_good_cells])), cex=0.8)
  grid(col='gray')
  dev.off()

  # global distrib of Pr(dup)
  png(sprintf("%s/emp_p_dup.png", sch_fig_dir), width=250, height=250)
  plot(density(unlist(dup)), main="Pr(dup)", xlab="prob per chrom", lwd=2)
  grid(col="lightgray")
  dev.off()

  # z_dup: z-score of p_dup vs expected by binomial multi-cell model

  chr_fe = table(gextract("redb.feGATC_gc", gintervals.all())$chrom)
  chr_fe = chr_fe[grep("chrY|chrM", names(chr_fe), invert=T)]

  if (!file.exists(sprintf("%s/cell_chr_cov.txt", sch_table_dir))) {
    sch_chrom_cov()
  }
  cov = read.table(sprintf("%s/cell_chr_cov.txt", sch_table_dir), header=T)
  cov = cov[rownames(dup), grep("chrY|chrM", colnames(cov), invert=T)]

  m_fe = matrix(chr_fe, nrow(cov), ncol(cov), byrow=T)

  # e_dup vs o_dup
  calc_e_dup = function(k, n) {
    (1 - ((k-1)/k)**n - (n/k)*((k-1)/k)**(n-1)) / (1 - ((k-1)/k)**n)
  }

  # compare expected trend to o_dup
  ex_chr = paste0("chr", sch_chroms[1])
  png(sprintf("%s/obs_vs_exp_dup.png", sch_fig_dir), width=350, height=350)
  plot(cov[,ex_chr], sapply(cov[,ex_chr], calc_e_dup, k=chr_fe[ex_chr]), type='p', pch=19, cex=0.2, main="chrom 1", xlab="#contacts", ylab="Pr(dup)", col='blue', ylim=range(dup[,ex_chr]))
  points(cov[grep('NXT', rownames(cov)),ex_chr], dup[grep('NXT', rownames(cov)),ex_chr], pch=19, cex=0.2, col='red')
  points(cov[grep('NST', rownames(cov)),ex_chr], dup[grep('NST', rownames(cov)),ex_chr], pch=19, cex=0.2, col='darkgreen')
  grid(col="lightgray")
  legend("bottomright", legend=c("Exp (binom)", "Obs (2i)", "Obs (serum)"), pch=19, col=c("blue", "red", "darkgreen"), cex=0.9, bty='n')
  dev.off()

  e_dup = calc_e_dup(k=m_fe, n=cov)
  #z_dup = (dup - e_dup) / sqrt(e_dup)

  # z-dup distrib
#  png(sprintf("%s/z_dup.png", sch_fig_dir), width=250, height=250)
  #plot(density(unlist(z_dup)), main="dup z-score", xlab="z-dup", lwd=2, xlim=c(min(z_dup), quantile(unlist(z_dup), 0.995)))
  #grid(col="lightgray")
  #dev.off()

  # dup vs tor marginal cov over early rep regions
  png(sprintf("%s/avg_dup_vs_tor_marg.png", sch_fig_dir), width=300, height=300)
  mu_d = rowMeans(dup)
  tor_s = rowSums(sch_tor_marg_n[rownames(dup), c('tor3', 'tor4')])
  plot(mu_d, tor_s, type='p', pch=19, cex=0.2, main="", xlab="mean dup per cell", ylab="early tor cov marg", xlim=c(min(dup), quantile(unlist(dup), 0.99)))
  legend("topright", legend=sprintf("rho=%.2f", cor(mu_d, tor_s)), pch=19)
  grid(col="lightgray")
  dev.off()

  # distrib of mean dup and early tor coverage by batch DNA sort groups
  if (!is.null(n1_b)) {
    n1_b = c(1845, 1994, 2005, 2014, 2027)
    n1_nms = intersect(rownames(dup), rownames(sch_batch)[is.element(sch_batch$batch_id, n1_b)])

    nons_b = c(1, 1268, 1527, 1632, 1654, 1749, 2434, 2435)
    nons_nms = intersect(rownames(dup), rownames(sch_batch)[is.element(sch_batch$batch_id, nons_b)])

    serum_b = c(3462, 3463, 3580, 3751, 3772, 3941)
    serum_nms = intersect(rownames(dup), rownames(sch_batch)[is.element(sch_batch$batch_id, serum_b)])

    png(sprintf("%s/dup_and_tor_cov_by_batch_grp.png", sch_fig_dir), width=350, height=600)
    layout(matrix(1:2, nrow=2, ncol=1))

    plot(density(rowMeans(dup[n1_nms,])), main="Pr(dup)", xlab="", col='blue', lwd=2)
    lines(density(rowMeans(dup[nons_nms,])), col='red', lwd=2)
    lines(density(rowMeans(dup[serum_nms,])), col='darkgreen', lwd=2)
    grid()
    legend("topright", legend=c("n1 sorted", "non sorted", "serum (<2n)"), lty=1, lwd=2, col=c('blue', 'red', 'darkgreen'), bty='n', cex=0.9)

    plot(density(rowSums(sch_tor_marg_n[n1_nms,c('tor3', 'tor4')])), main="Early/Total coverage", xlab="", col='blue', lwd=2)
    lines(density(rowSums(sch_tor_marg_n[nons_nms,c('tor3', 'tor4')])), col='red', lwd=2)
    lines(density(rowSums(sch_tor_marg_n[serum_nms,c('tor3', 'tor4')])), col='darkgreen', lwd=2)
    grid()
    legend("topright", legend=c("n1 sorted", "non sorted", "serum (<2n)"), lty=1, lwd=2, col=c('blue', 'red', 'darkgreen'), bty='n', cex=0.9)

    dev.off()
  }


  # dup chr1 vs all
  png(sprintf("%s/dup_%s_vs_all.png", sch_fig_dir, ex_chr), width=1200, height=900)
  par(mar=c(7,6,6,0))
  layout(matrix(seq(1,20,1),nrow=4,ncol=5))
  par(mar=c(1,1,1,1))

  cols = colnames(dup)
  for(col in cols[-1]) {
    plot(dup[,cols[1]], dup[, col], pch=19, xlim=c(0,0.04), ylim=c(0,0.04), cex=0.5)
    legend("topleft", legend=sprintf("%s (rho=%.2f)", col, cor(dup[,cols[1]], dup[,col])), bty='n', pch=19, cex=1.5)
    grid()
  }
  dev.off()

  # z-dup by batch
  good_cells_batch = sch_batch[intersect(rownames(dup), sch_good_cells), ]
  good_cells_batch = good_cells_batch[order(good_cells_batch$batch), ]

  batch_numbers = unique(good_cells_batch$batch_id)
  ind = c(1, 1+cumsum(rle(good_cells_batch$batch)$lengths))
  ind = ind[-length(ind)]

  dupg = dup[rownames(good_cells_batch),]

  png(sprintf("%s/dup_by_batch.png", sch_fig_dir), width=1500, height=600)
  par(mar=c(7,6,6,0))
  plot(rowMeans(dupg), pch=21,
       bg=sch_batch_colors[good_cells_batch$batch],
       col="black",
       xaxt="n",
       xlab="",
       ylab="",
       cex.axis=2)
  #title(xlab="Batch number", cex.lab=2, line=5)
  title(ylab="Mean Fend-Dup", cex.lab=2, line=4)
  abline(v=ind, col="gray")
  med = median(rowMeans(dup))
  abline(h=med)
  axis(side=1, at=ind[c(TRUE, FALSE)], labels=batch_numbers[c(TRUE, FALSE)], las=2, cex.axis=1)
  axis(side=3, at=ind[c(FALSE, TRUE)], labels=batch_numbers[c(FALSE, TRUE)], las=2, cex.axis=1)
  #points(apply(dupg, 1, min), pch=21, bg="gray", col="black", cex=0.4)
  #points(apply(dupg, 1, max), pch=21, bg="gray", col="black", cex=0.4)
  dev.off()

  # cluster z-dup matrix by chromosomes and cells
  png(sprintf("%s/dup_clust.png", sch_fig_dir), width=1200, height=500)
  layout(matrix(2:1, nrow=2, ncol=1), height=c(1,9))
  par(mar=c(1,1,1,1))

  hc = hclust(dist(dup), "ward.D2")
  hc_chr = hclust(dist(t(dup)), "ward.D2")
  breaks = unique(c(seq(min(dup), med, length=51), seq(med, max(dup), length=51)))
  image(as.matrix(dup[hc$order,hc_chr$order]), col=colorRampPalette(c("blue", "black", "red"))(100), breaks=breaks, xaxt="n", yaxt='n', xlab="good cells", ylab="")
  axis(2, at=seq(0, 1, length=ncol(dup)), labels=gsub("chr", "", colnames(dup)[hc_chr$order]), tick=F, las=1, cex=0.8)

  mdup = rowMeans(dup[hc$order,])
  image(as.matrix(rowMeans(dup[hc$order,]) > max_dup), col=colorRampPalette(c("white", "red"))(100), xlab="", ylab="", xaxt='n', yaxt='n' )

  dev.off()

}


##################
calc_ab_enrich <- function(ab, name, do_plots=F) {
  colnames(ab) = toupper(gsub(paste(name, ".", sep=""), "", colnames(ab)))

  pA = (2*ab$AA + ab$AB) / (2 * rowSums(ab))
  pB = (2*ab$BB + ab$AB) / (2 * rowSums(ab))

  ab$enAA = log2( ab$AA / (pA**2 * rowSums(ab)))
  ab$deplAB = -log2( ab$AB / (2*pA*pB * rowSums(ab)))
  ab$enBB = log2( ab$BB / (pB**2 * rowSums(ab)))

  if (do_plots) {
    png(sprintf("%s/ab_%s_enrich_distr.png", sch_fig_dir, name), 300, 300)
    plot(density(ab$enAA), xlim=range(ab[,c('enAA', 'enAB', 'enBB')]), col='blue', lwd=2, main=sprintf("%s A,B", name), xlab="log2(o/e)")
    lines(density(ab$enAB), col='orange', lwd=2)
    lines(density(ab$enBB), col='darkred', lwd=2)
    grid(col="gray")
    legend("topleft", legend=c("AA", "AB", "BB"), lwd=2, lty=1, col=c("blue", "orange", "darkred"), cex=0.8, bty='n')
    dev.off()
  }
  ab
}



#===================================
#
sch_plot_mean_mat <- function(base_dir, clusts, no_dig_dist=1000, w_per_100m = 200, dens_white_quantile=0.1, scope = gintervals.2d.all(), force_size = -1, bwx=1e5, bwy=bwx, fn_pref="", ngrid=1000, colspec = c("white", "white", "white", "white", "white", "white", "white", "white", "white", "white", "gray", "gray", "lightblue", "blue", "darkblue", "orange", "yellow", "black"), glob_dens_lim=NULL)
{

    shades = colorRampPalette(colspec)
    cl_idxs = list()
    for (cl in 1:max(clusts$x)) {
      cl_idxs[[cl]] = which(clusts$x==cl)
    }

    names = paste0("C", 1:max(clusts$x), " ", table(clusts$x))

    commands = paste0(".sch_plot_group_mean_mat(base_dir, rownames(clusts)[cl_idxs[[", 1:max(clusts$x), "]]], no_dig_dist=no_dig_dist, w_per_100m=w_per_100m, dens_white_quantile=dens_white_quantile, scope=scope, force_size=force_size, bwx=bwx, bwy=bwy, fn_pref=\"", paste0(fn_pref, 1:max(clusts$x)), "\", ngrid=ngrid, glob_dens_lim=glob_dens_lim, name=\"", names, "\", shades)", collapse=",")

    res = eval(parse(text=paste("gcluster.run(",commands,")")))
  }

.sch_plot_group_mean_mat <- function(base_dir, nms, no_dig_dist=1000, w_per_100m = 200, dens_white_quantile=0.1, scope = gintervals.2d.all(), force_size = -1, bwx=1e5, bwy=bwx, fn_pref="", ngrid=1000, shades=colorRampPalette(c("white", "white", "white", "white", "white", "white", "white", "white", "white", "white", "gray", "gray", "lightblue", "blue", "darkblue", "orange", "yellow", "black")), name="", type="kde", add_axis=F, point_cex=0.5, rotate=F, glob_dens_lim=NULL, ablines_x=NULL, ablines_y=NULL, ablines_lwd=1)
{
  library(KernSmooth)

  w_per_bp = w_per_100m/1e+8
  ag1 = gintervals.all()

  mat = gextract(nms[1], scope, colnames=c("obs"))
  for (nm in nms[-1]) {
    #message("will extract ", nm)
    cur_mat = gextract(nm, scope, colnames=c("obs"))
    cur_mat = cur_mat[abs(cur_mat$start1 - cur_mat$start2) > no_dig_dist,]
    mat = rbind(mat,cur_mat)
  }
  message("mat has ", nrow(mat), " entries")

  for (chr_nm in unique(scope$chrom1)) {
    if (chr_nm == "chrM" || chr_nm == "chrY") {
      next
    }
    if (force_size == -1) {
      sz = round(ag1$end[ag1$chrom==chr_nm]*w_per_bp)
    } else {
      sz = force_size
    }
    if (nrow(scope) == 1) {
      mat_c = mat
    } else {
      mat_c = mat[mat$chrom1 == chr_nm & mat$chrom2 == chr_nm,]
    }
    if (nrow(mat_c) < 5) {
      message("no data for chr ", chr_nm)
      next
    }
    system(paste("mkdir -p", base_dir))

    fig_png(sprintf("%s/%s_%s.png", base_dir, fn_pref, chr_nm), width=sz, height=sz)
    par(mar=c(1,1,3,1))
    if (rotate) {
      stopifnot(type == "points", scope$start1 == scope$start2, scope$end1 == scope$end2)

      xlim = c(scope$start1, scope$end1)
      ylim = c(0, scope$end1 - scope$start1)/2

      mat = mat[ mat$start2 > mat$start1,]
      y = (mat$start1 + mat$start2)/2
      x = (mat$start2 - mat$start1)/2
    }
    else {
      xlim=c(scope$start1, scope$end1)
      ylim=c(scope$start2, scope$end2)
      x = mat$start1
      y = mat$start2
    }

    if (type == "kde") {
      d2d = bkde2D(cbind(x, y),
        bandwidth=c(bwx,bwy),
        gridsize=rep(ngrid, 2))
      dens = d2d$fhat
      lowdens = quantile(dens[dens>0], dens_white_quantile, na.rm=T)
      highdens = quantile(dens[dens>0], 0.999,na.rm=T)
      if (!is.null(glob_dens_lim)) {
        lowdens = glob_dens_lim[1]
        highdens = glob_dens_lim[2]
      }
      message(sprintf("dens lim %e and %e", lowdens, highdens))
      dens[dens<lowdens] = lowdens
      dens[dens>highdens] = highdens

#      image(d2d$x1, d2d$x2, t(log10(dens+lowdens)), col=shades(1000),
      image(d2d$x1, d2d$x2, t(log10(dens)), col=shades(1000),
            zlim=c(log10(lowdens*2), log10(highdens+lowdens)), useRaster=TRUE,
            xaxt="n", yaxt="n", xlab="", ylab="", main="", xlim=xlim, ylim=ylim , asp=1)
    }
    else if (type == "points") {
      if (rotate) {
        x = (mat$start1 + mat$start2)/2
        y = (mat$start2 - mat$start1)/2

        x = x[y > 0]
        y = y[y > 0]
      }
      else {
        x = mat$start1
        y = mat$start2
      }

      plot(x, y, type='p', pch=19, cex=point_cex, xlab="", ylab="", xaxt="n", yaxt="n", xlim=xlim, ylim=ylim, asp=1)
    }

    title(paste(chr_nm, name, nrow(mat)))

    if (add_axis) {
      axis(1)
      axis(2)
    }


    if (!is.null(ablines_x)) {
      abline(a=0, b=1, lwd=ablines_lwd)
      abline(v=ablines_x, lwd=ablines_lwd)
    }
    if (!is.null(ablines_y)) {
      abline(a=0, b=1, lwd=ablines_lwd)
      abline(h=abines_y, lwd=ablines_lwd)
    }
    dev.off()
  }

  plot_legend(sprintf("%s/%s_%s_lim%g_to_%g_leg.png", base_dir, fn_pref, chr_nm, log10(glob_dens_lim[1]), log10(glob_dens_lim[2])), log10(glob_dens_lim), shades(1000), "density")

}

##############################################
# Pooling contacts from single cells
# Binarizing contact count before summing, so "f1 f2 x" means contact f1-f2 appeared in x cells
sch_create_pooled_track <- function(cells_dir, pool_track_nm=pool_tn, adj_base_dir="/net/mraid14/export/data/users/lubling/datasets/scell", redb=sch_redb_dir)
{
  if (is.null(adj_base_dir)) {
    adj_files = paste(cells_dir, "adj", sep="/")
  }
  else {
    adj_files = paste(adj_base_dir, cells_dir, "adj", sep="/")
  }

  gtrack.2d.import_contacts(pool_track_nm, description="pooled scell adj", contacts=adj_files, fends=paste(redb, "/GATC.fends", sep=""))
}

#==============================================================================
# "sch_create_pooled_track_wrapper" creates pooled tracks
#==============================================================================

sch_create_pooled_track_wrapper <- function(nms, pool_nm)
{
  if (haploid) {
    all_nxt_dirs = paste(sch_base_dir, "2i", gsub("_", ".", gsub("scell.nextera.", "", grep("NXT", nms, value=T)), fixed=T), sep="/")
    all_nst_dirs = paste(sch_base_dir, "serum", gsub("_", ".", gsub("scell.nextera.", "", grep("NST", nms, value=T)), fixed=T), sep="/")
    
    if (pool_tn == "scell.nextera.pool_good_hap_2i_serum_es") {
      all_dirs = c(all_nxt_dirs, all_nst_dirs)
    }
    else if (pool_tn == "scell.nextera.pool_good_hap_2i_es") {
      all_dirs = all_nxt_dirs
    }
    else if (pool_tn == "scell.nextera.pool_good_hap_serum_es") {
      all_dirs = all_nst_dirs
    }
  }
  else {
    nms_regexp = switch(pool_tn,
      scell.nextera.pool_good_hyb_2i_sort_es="1CDX|1CDU",
      scell.nextera.pool_good_hyb_serum_es="1CDS",
      scell.nextera.pool_good_hyb_2i_all_es="1CDX|1CDU|1CDES",
      scell.nextera.pool_good_hyb_2i_no_sort_es="1CDES")

    all_names = gsub("_", ".", gsub("scell.nextera.", "", grep(nms_regexp, nms, value=T, perl=T)), fixed=T)
    all_names = gsub("1CDES.", "1CDES_", all_names, fixed=T)
    
    all_dirs = unlist(sapply(sch_base_dir, function(x) { paste(x, all_names, sep="/") }))
    all_dirs = all_dirs[file.exists(all_dirs)]

  }

  message(sprintf("======== creating pool track %s... from %d files" ,pool_nm, length(all_dirs)))
  sch_create_pooled_track(all_dirs, pool_track_nm=pool_nm, adj_base_dir="")
}


################
#
sch_gen_group_pools <- function(g_nms=c("post_m", "g1", "early_s", "mid_s_g2", "pre_m"))
{
  nms = split(rownames(sch_decay_metrics), sch_decay_metrics$group)
  g_tns = paste(pool_tn, "group", seq_along(g_nms), g_nms, sep="_")
  
  for (i in 1:length(g_nms)) {
    if (gtrack.exists(g_tns[i])) {
      gtrack.rm(g_tns[i], force=T)
    }
    message("writing track ", g_tns[i])
    sch_create_pooled_track_wrapper(nms[[i]], g_tns[i])
  }

}

#########
.sch_get_pool_domains <- function(dom_cutoff_q=ins_dom_thresh, res=1e3, add_epi_tracks=T, discard_below=1000, ins_track=ins_tn, raw_tn=pool_tn)
{
  if (!gtrack.exists(ins_track)) {
    gtrack.2d.gen_insu_track(raw_tn, ins_scale, res=res, new_track=ins_track, min_diag_d=discard_below)
  }
  cd = gtrack.2d.get_insu_doms(ins_track, gquantiles(ins_track, dom_cutoff_q), iterator=res)
  cd$len = cd$end - cd$start
  cd = cd[cd$len >= domain_size_range[1] & cd$len <= domain_size_range[2], ]

  na_d = gscreen(sprintf("is.na(%s)", ins_track), iterator=res)
  d = gintervals.neighbors(cd, na_d[na_d$end - na_d$start >= intra_domain_na_segment_len_to_filter, ], maxneighbors=1, maxdist=0, na.if.notfound=T)
  d = d[is.na(d[,5]), 1:4]

  d$chrom = factor(d$chrom, levels=paste0("chr", sch_chrom_levels))
  d = d[order(as.numeric(d$chrom)),]

  d$id = paste(d$chrom, d$start, sep="_")
  rownames(d) = paste(d$chrom, d$start, sep="_")

  if (add_epi_tracks) {
    tor = gextract(sch_tor_track, intervals=d, iterator=d, colnames="tor")
    tor[order(tor$intervalID),"id"] = d$id
    d[tor$id,"tor"] = tor$tor
  }

  d
}

#########
.sch_get_pool_tad_borders <- function(iterator=1e3, input_track=pool_tn, ins_track=ins_tn, scale=ins_scale, min_diag=sch_remove_near_dist, filter_by_cov_per_kb=50, filter_iter=0, scope=NULL)
{
  if (!gtrack.exists(ins_track)) {
    gtrack.2d.gen_insu_track(input_track, scale, res=iterator, min_diag_d=min_diag, new_track=ins_track)
  }
  message("got ins track")
  ins = gscreen(sprintf("!is.na(%s) & %s <= %f", ins_track, ins_track, gquantiles(ins_track, ins_dom_thresh)), iterator=iterator)
  message("found insulating regions")
  chroms = gintervals.all()
  rownames(chroms) = as.character(chroms$chrom)
  ins$max = chroms[ as.character(ins$chrom), 'end']
  ins$start = pmax(ins$start - extend_border_region, 0)
  ins$end   = pmin(ins$end   + extend_border_region, ins$max)
  insc = gintervals.canonic(ins, unify=T)

  gvtrack.create("ins_min", ins_track, "min")
  bords = gextract("ins_min", intervals=insc, iterator=iterator)
  message("extracted min from region")
  bp = do.call(rbind, lapply(split(bords, bords$intervalID), function(x) { x[which.min(x$ins)[1],]}))

  bp = intervals.centers(bp)

  if (filter_iter > 0) {
    if (is.null(scope)) {
      scope = .get_well_covered_scope(filter_iter, pool_tn, filter_by_cov_per_kb)
    }
    bp = gintervals.intersect(bp, scope)
  }

  bp
}

##############################################
# Cluster domains by inter-domain profile
cluster_pool_hic_domains_by_inter_domain_contacts <- function(cluster_method='kmeans', k=2, cluster_by='trans', scale_by='area', delta=1e-10, trim_below=-5, breaks_mode="around.median", q_disp=0.01, heatmap_cols=c("white", "darkred"), fn_pref=NULL, do_plot=T, vfunc="area", d=NULL, tn=pool_tn, ins_track=ins_tn, add_epi_data=T)
{
  if (is.null(d)) {
    d = .sch_get_pool_domains(raw_tn=tn, ins_track=ins_track)
  }

  if (add_epi_data) {
    if (gtrack.exists(active_tn)) {
      gvtrack.create("k4me3", active_tn, "global.percentile.max")
      k4peaks = gextract("ifelse(k4me3 >= 0.99, 1, 0)", intervals=d, iterator=50, colnames="k4me3")
      d$k4me3 = tapply(k4peaks$k4me3, k4peaks$intervalID, mean)
    }
    if (gtrack.exists(lamin_tn)) {
      d$laminB1 = gextract(lamin_tn, intervals=d, iterator=d, colnames="lamin")$lamin
    }
    if (gtrack.exists("Encode.esb4.h3k27me3.rep1") & gtrack.exists("Encode.esb4.h3k27me3.rep2")) {
      d$k27me3 = gextract("(Encode.esb4.h3k27me3.rep1 + Encode.esb4.h3k27me3.rep2)/2", intervals=d, iterator=d, colnames='x')$x
    }
  }
  
  cat("Creating intervals ...\n")
  d_inds = as.data.frame(expand.grid(1:nrow(d), 1:nrow(d)))
  d_ints = data.frame(chrom1=d[d_inds[,1],'chrom'], start1=d[d_inds[,1],'start'], end1=d[d_inds[,1],'end'], chrom2=d[d_inds[,2],'chrom'], start2=d[d_inds[,2],'start'], end2=d[d_inds[,2],'end'])
  ints = gintervals.canonic(d_ints)

  ps_tn = paste(tn, vfunc, sep="_")
  gvtrack.create(ps_tn, tn, vfunc)

  cat(sprintf("Extracting from %s ...\n", ps_tn))
  d_cov = gextract(sprintf("ifelse(is.na(%s), 0, %s)", ps_tn, ps_tn) , intervals=ints, iterator=ints, colnames="cov")

  if (cluster_by == 'trans') {
    d_cov[ as.character(d_cov$chrom1) == as.character(d_cov$chrom2), 'cov'] = NA
  }
  m = matrix(nrow=nrow(d), ncol=nrow(d), data=0)
  d_reord_inds = d_inds[attr(ints, 'mapping'),]

  m[cbind(d_reord_inds[,1], d_reord_inds[,2])] = d_cov$cov
  rownames(m)[d_reord_inds[,1]] = paste(d_cov$chrom1, d_cov$start1, sep="_")

  if (scale_by == 'simple-mean') {
    m = log2((m+delta) / mean(m+delta, na.rm=T))
  }
  else if (scale_by == 'area') {
    s = matrix(nrow=nrow(m), ncol=ncol(m), data=0)
    s[ cbind(d_inds[,1], d_inds[,2]) ] = d[ d_inds[,1], 'len'] * d[ d_inds[,2], 'len']
    s = (s * sum(m+delta, na.rm=T)) / sum(s[!is.na(m)])
    m = log2( (m+delta) / s)
  }
  else if (scale_by == 'marginal') {
    marg = rowSums(m+1, na.rm=T) / sum(m+1, na.rm=T)
    s = matrix(nrow=nrow(m), ncol=ncol(m), data=0)
    s[ cbind(d_inds[,1], d_inds[,2]) ] = marg[ d_inds[,1]] * marg[ d_inds[,2]]
    s = s * sum(m+delta, na.rm=T)
    m = log2( (m+delta) / s)
  }
  if (is.null(fn_pref)) {
    fn_pref = sch_ab_cluster_fn_pref
    # fn_pref = sprintf("%s/pool_domains_cluster_%s_%s_%s_%s_k%d_d%f_breaks.%s_q.disp%f_%dcols_min%d", sch_table_dir, cluster_by, scale_by, vfunc, cluster_method, k, signif(log10(delta), 2), breaks_mode, q_disp, length(heatmap_cols), trim_below)
  }

  if (do_plot) {
    png(sprintf("%s_hist.png", fn_pref), width=420, height=420)
    hist(m, 150, main=scale_by, xlab="contact enrichment (log2)")
    grid(col="gray")
    abline(v=trim_below, col='red')
    dev.off()
  }

  if (!is.na(trim_below)) {
    m[ m < trim_below ] = trim_below
  }

  cat(sprintf("%s clustering...\n", cluster_method))
  cx = hc = km = NULL
  if (cluster_method == 'hclust') {
    hc = hclust(dist(m))
    cx = m[hc$order, hc$order]
  }
  else if (cluster_method == 'kmeans') {
    km = TGLKMeans_wrapper(m, sprintf("%s.tab", fn_pref), k)
    hc = hclust(dist(km$centers))

    new_order = c()
    for (i in hc$order) {
      new_order = c(new_order, which(km$cluster == i))
    }
    cx = m[new_order, new_order]

    d$cluster = km$cluster
  }

  if (!is.na(q_disp)) {
    quants = quantile(unlist(cx), c(q_disp, 1 - q_disp), na.rm=T)
    cx [ cx < quants[1] ] = quants[1]
    cx [ cx > quants[2] ] = quants[2]
  }

  if (do_plot) {
    my_palette = colorRampPalette(heatmap_cols)(n = (length(heatmap_cols)-1)*100 - 1)

    if (length(heatmap_cols) == 3) {
      col_breaks = switch(breaks_mode,
        around_zero=c(seq(min(cx, na.rm=T), 0, length=100), seq(0, max(cx, na.rm=T), length=100)),
        around_median=c(seq(min(cx, na.rm=T), median(cx, na.rm=T), length=100), seq(median(cx, na.rm=T), max(cx, na.rm=T), length=100)))
    }
    else if (length(heatmap_cols) == 2) {
      col_breaks = c(seq(min(cx, na.rm=T), max(cx, na.rm=T), length=100))
    }

    png(sprintf("%s_legend.png", fn_pref), width=200, height=400)
    plot.new()
    plot.window(xlim=c(0,5), ylim=c(1,100))
    rect(xleft=1, xright=2, ybottom=1:99, ytop=2:100, col=my_palette, border=NA)
    rect(xleft=1, xright=2, ybottom=1, ytop=100, border="black")
    segments(x0=2, x1=2.1, y0=c(1,25, 50, 75, 100), y1=c(1,25, 50, 75, 100))
    text(x=2.2, y=c(1,50,100), labels=signif(col_breaks[c(1,50,100)], 2), pos=4)
    dev.off()

    png(sprintf("%s.png", fn_pref), width=nrow(cx)+50, height=ncol(cx)+50)
    image(as.matrix(cx), col=my_palette, breaks=col_breaks, main="", xlab="", ylab="", xaxt='n', yaxt='n')
    
    if (cluster_method == 'kmeans') {
      clus_lines = cumsum(km$size[hc$order]) / ncol(cx)
      clus_mid = (clus_lines + c(0, clus_lines[-length(clus_lines)]))/2
      
      axis(2, at=clus_mid, labels=paste(hc$order, " #", km$size[hc$order], sep=""), tick=F, las=1)
      abline(h=clus_lines, lw=6)
      abline(v=clus_lines, lw=6)
    }
    
    dev.off()
    
    png(sprintf("%s_pool_domains_cluster_genomic_%s.png", fn_pref, scale_by), width=nrow(cx)+50, height=ncol(cx)+50)
    
    image(as.matrix(m), col=my_palette, breaks=col_breaks, main="", xlab="", ylab="", xaxt='n', yaxt='n')
    
    dev.off()
  }

  list(m=m, cx=cx, km=km, hc=hc, d=d)
}


############
#
sch_plot_genome_wide <- function(nm, ofn, chroms=sch_chroms, trans_cols=c('white', '#ffffcc','#ffeda0', '#fed976', '#feb24c'), cis_add_cols=c('#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026', 'purple', 'black', 'white'), binsize=4e6, min_contacts=0, zlim=NULL, reg=1e-5, vfunc="area", legend_ofn=NULL)
{
  options(error=dump.frames)
  chrs = gintervals.all()
  rownames(chrs) = chrs$chrom
  chrs = chrs[paste("chr", chroms, sep=""), ]

  chrs$glob_start = cumsum(c(0, chrs$end[-nrow(chrs)]))
  chrs$glob_end = cumsum(chrs$end)

  chrs2d = gintervals.2d.all()
  chrs2d = chrs2d[is.element(chrs2d$chrom1, chrs$chrom) & is.element(chrs2d$chrom2, chrs$chrom), ]

  iter1d = giterator.intervals(intervals=gintervals(chroms), iterator=binsize)
  iind = expand.grid(1:nrow(iter1d), 1:nrow(iter1d))
  i1 = iind[,1]
  i2 = iind[,2]
  iter2d = gintervals.2d(iter1d$chrom[i1], iter1d$start[i1], iter1d$end[i1], iter1d$chrom[i2], iter1d$start[i2], iter1d$end[i2])

  vtn = paste(nm, vfunc, sep="")
  gvtrack.create(vtn, nm, vfunc)

  # get contacts
  mat = gextract(sprintf("ifelse(!is.nan(%s),%s,0)", vtn, vtn), intervals=chrs2d, iterator=iter2d, colnames='count')

  mat = merge(mat, chrs[,c('chrom', 'glob_start')], by.x="chrom1", by.y="chrom")
  mat = merge(mat, chrs[,c('chrom', 'glob_start')], by.x="chrom2", by.y="chrom")

  mat$start1 = mat$start1 + mat$glob_start.x
  mat$start2 = mat$start2 + mat$glob_start.y

  mat[mat$count < min_contacts, 'count'] = NA
  mat$count = log2(mat$count + reg)
  cis = as.character(mat$chrom1) == as.character(mat$chrom2)

  if (is.null(zlim)) {
    zlim = c(0, quantile(mat[cis, 'count'], 0.99, na.rm=T))
  }

  # plot cmap
  if (!is.null(ofn)) {
    fig_png(ofn, width=nrow(iter1d) + 80, height=nrow(iter1d) + 100)
  }
  par(mar=c(3,3,3,1))
  message("plotting...")
  starts = sort(unique(mat$start1))
  m = matrix(NA, length(starts), length(starts))
  rownames(m) = colnames(m) = starts
  m[ cbind(as.character(mat$start1), as.character(mat$start2)) ] = mat$count

  mid_val = max(0, min(quantile(mat[!cis, 'count'], 0.99, na.rm=T), quantile(mat[cis, 'count'], 0.01, na.rm=T)))

  #print(c(zlim[1], mid_val, zlim[2]))
  cols = diverge_color(c(trans_cols, cis_add_cols), c(seq(zlim[1], mid_val, len=length(trans_cols)), seq(mid_val, zlim[2], len=length(cis_add_cols))), 1000)
  image(starts + binsize/2, starts + binsize/2, pmin(pmax(m, zlim[1]), zlim[2]), col=cols, zlim=zlim, xaxt="n", yaxt="n", xlab="", ylab="", main=nm)

  clines = chrs$glob_end
  abline(h=clines, col='black', lwd=1)
  abline(v=clines, col='black', lwd=1)

  clabpos = (chrs$glob_start + chrs$glob_end) / 2
  axis(1, at=clabpos, labels=gsub("chr", "", chrs$chrom), tick=F, cex=2)
  axis(2, at=clabpos, labels=gsub("chr", "", chrs$chrom), tick=F, cex=2)

  if (!is.null(ofn)) {
    dev.off()
  }

  if (!is.null(legend_ofn)) {
    plot_legend(legend_ofn , zlim, cols, "log2(count)")
  }
  mat
}

#############
#
plot_trans_helper <- function(nm, chroms, cex=0.3, alpha=0.5, bp_per_px=5e5, odir=NULL, grid_every_bp=20e6)
{

  chrs = gintervals.all()
  rownames(chrs) = chrs$chrom
  chrs = chrs[paste("chr", chroms, sep=""), ]

  chrs$glob_start = cumsum(c(0, chrs$end[-nrow(chrs)]))
  chrs$glob_end = cumsum(chrs$end)

  chrs2d = gintervals.2d.all()
  chrs2d = chrs2d[is.element(chrs2d$chrom1, chrs$chrom) & is.element(chrs2d$chrom2, chrs$chrom), ]

  mat = gextract(nm, intervals=chrs2d, colnames='count')

  mat = merge(mat, chrs[,c('chrom', 'glob_start')], by.x="chrom1", by.y="chrom")
  mat = merge(mat, chrs[,c('chrom', 'glob_start')], by.x="chrom2", by.y="chrom")

  mat$start1 = mat$start1 + mat$glob_start.x
  mat$start2 = mat$start2 + mat$glob_start.y

  if (!is.null(odir)) {
    size = max(chrs$glob_end) / bp_per_px
    fig_png(sprintf("%s/%s_chrs_%s.png", odir, nm, paste0(chroms, collapse="_")), size, size, pointsize=pointsize, res=fres)
  }
  par(mar=c(2,2,2,0.5))


  plot(mat$start1, mat$start2, pch=19, cex=cex, col=addalpha('blue', alpha), main=sprintf("%s (g%s, %d)", sub("scell.nextera.", "", nm), sch_decay_metrics[nm, 'group'], sch_decay_metrics[nm, 'ord']), xlim=c(0, max(chrs$glob_end)), ylim=c(0, max(chrs$glob_end)), xaxt='n', yaxt='n', xlab="", ylab="")

  clines = c(0, chrs$glob_end)
  abline(h=clines, col='black', lwd=2)
  abline(v=clines, col='black', lwd=2)

  clabpos = (chrs$glob_start + chrs$glob_end) / 2
  axis(1, at=clabpos, labels=gsub("chr", "", chrs$chrom), tick=F, cex=2)
  axis(2, at=clabpos, labels=gsub("chr", "", chrs$chrom), tick=F, cex=2)

  if (grid_every_bp > 0) {
    grids = do.call("c", sapply(rownames(chrs), function(chr) { seq(chrs[chr, 'glob_start'], chrs[chr, 'glob_end'], by=grid_every_bp) }))

    abline(h=grids, col='black', lwd=0.5, lty=2)
    abline(v=grids, col='black', lwd=0.5, lty=2)
  }

  if (!is.null(odir)) {
    dev.off()
  }

}


###########
count_dup_trans_contacts_across_all_cells <- function()
{
  message("count_dup_trans_contacts_across_all_cells")
  # heavy dups are assembly problems (e.g. chr11 genome-wide line)
  chroms = gintervals.2d.all()
  chroms = chroms[as.character(chroms$chrom1) < as.character(chroms$chrom2), ]
  chroms = chroms[!is.element(chroms$chrom1, c('chrM', 'chrY')) & !is.element(chroms$chrom2, c('chrM', 'chrY')), ]
  commands = paste0("count_dup_trans_contacts_in_spec_chroms_across_all_cells(chr1=\"", chroms$chrom1, "\", chr2=\"", chroms$chrom2, "\")", collapse=",")
  res = eval(parse(text=paste("gcluster.run(",commands,")")))
  maxc = max(unlist(lapply(res, function(x) max(as.numeric(names(x$retv$count))))))

  stats = rep(0, maxc)
  mul = NULL
  for (r in res) {
    counts = r$retv$count
    stats[as.numeric(names(counts))] = stats[as.numeric(names(counts))] + counts
    mul = rbind(mul, r$retv$mul)
  }

  sch_trans_switchers <<- data.frame(n=1:length(stats), count=stats)

  write.table(sch_trans_switchers, sprintf("%s/trans_dup_contacts.txt", sch_table_dir), quote=F, sep="\t")
  mul
}

count_dup_trans_contacts_in_spec_chroms_across_all_cells <- function(chr1="chr15", chr2="chr16")
{
  nms = gtrack.ls(sch_track_base)
  ints = gintervals.2d(chr1, 0, -1, chr2, 0, -1)
  tab = table(do.call("c", lapply(nms, function(nm) {  x = gextract(nm, intervals=ints); if (!is.null(x)) { message(paste(nm, nrow(x))); paste(x$start1, x$start2) }})))

  if (sum(tab > 1) == 0) {
    df_mul = NULL
  }
  else {
    mul = tab[ tab > 1]
    s1 = as.numeric(gsub(" .*", "", names(mul)))
    s2 = as.numeric(gsub(".* ", "", names(mul)))
    df_mul = data.frame(chrom1=chr1, start1=s1, end1=s1+1, chrom2=chr2, start2=s2, end2=s2+1, count=mul)
  }
  list(count=table(tab), mul=df_mul)
}


########
compute_reads_per_contact <- function()
{
  message("compute_reads_per_contact")
  nms = gtrack.ls(sch_track_base)
  batches = split(rownames(sch_batch[nms, ]), sch_batch$batch)

  commands = paste0(".sch_track_reads_per_contact(nm=\"", nms, "\")", collapse=",")
  res = eval(parse(text=paste("gcluster.run(",commands,")")))
  stopifnot(sum(unlist(lapply(res, function(r) { typeof(r$retv) == 'list'} ))) == length(res))

  maxc = max(unlist(lapply(res, function(x) max(as.numeric(c(names(x$retv$cis), names(x$retv$trans)))))))

  cis_s = matrix(0, length(nms), maxc)
  rownames(cis_s) = nms
  trans_s = cis_s

  for (i in seq_along(nms)) {
    cis = res[[i]]$retv$cis
    trans = res[[i]]$retv$trans
    cis_s[nms[i], as.numeric(names(cis))] = cis
    trans_s[nms[i], as.numeric(names(trans))] = trans
  }

  sch_trans_contact_mul <<- trans_s
  sch_cis_contact_mul <<- cis_s

  write.table(sch_trans_contact_mul, sprintf("%s/sch_trans_contact_mul.txt", sch_table_dir), quote=F, sep="\t")
  write.table(sch_cis_contact_mul, sprintf("%s/sch_cis_contact_mul.txt", sch_table_dir), quote=F, sep="\t")

}

.sch_track_reads_per_contact <- function(nm)
{
  ints = gintervals.2d.all()
  cis = ints[ints$chrom1 == ints$chrom2, ]
  trans = ints[ints$chrom1 != ints$chrom2, ]

  xcis = gextract(nm, intervals=cis, colnames='cov')
  xtrans = gextract(nm, intervals=trans, colnames='cov')
  list(mu_cis=mean(xcis$cov), mu_trans=mean(xtrans$cov), cis=table(xcis$cov), trans=table(xtrans$cov))
}

########
compute_chrom_pairs_total_contacts <- function()
{
  message("compute_chrom_pairs_total_contacts")
  nms = gtrack.ls(sch_track_base)

  commands = paste0(".sch_track_chrom_pairs_contacts(nm=\"", nms, "\")", collapse=",")
  res = eval(parse(text=paste("gcluster.run(",commands,")")))
  stopifnot(sum(unlist(lapply(res, function(r) { typeof(r$retv) == 'double'} ))) == length(res))


  sch_trans_pairs_contact_count <<- do.call("rbind", lapply(res, function(r) r$retv))
  rownames(sch_trans_pairs_contact_count) = nms

  write.table(sch_trans_pairs_contact_count, sprintf("%s/sch_trans_pairs_contact_count_mul.txt", sch_table_dir), quote=F, sep="\t")

}

.sch_track_chrom_pairs_contacts <- function(nm)
{
  ints = gintervals.2d.all()
  trans_ints = ints[ints$chrom1 != ints$chrom2, ]
  gvtrack.create("vtn", nm, "area")

  xtrans = gextract("ifelse(is.na(vtn), 0, vtn)", intervals=trans_ints, iterator=trans_ints, colnames='cov')
  count = xtrans$cov
  names(count) = paste(xtrans$chrom1, xtrans$chrom2)
  count
}


####
fends_pair_orientation <- function(mindist=100, maxdist=1e4, chroms=sch_chroms, log_step=0.125)
{
  fends = gextract("redb.feGATC_flen", intervals=gintervals(chroms))

  bins = seq(log_step * floor(log2(mindist) / log_step), log_step * floor(log2(maxdist) / log_step), by=log_step)
  stats = matrix(0, 4, length(bins))
  rownames(stats) = c('f00', 'f01', 'f10', 'f11')
  colnames(stats) = paste0('b', bins)

  nms = gtrack.ls(sch_track_base)
  
  commands = paste0(".sch_fends_pair_orientation(nm=\"", nms, "\", fends=fends, mindist=mindist, maxdist=maxdist, chroms=chroms, log_step=log_step)", collapse=",")
  res = eval(parse(text=paste("gcluster.run(",commands,", max.jobs=300)")))

  for (i in seq_along(res)) {
    rv = res[[i]]$retv
    if (typeof(rv) == 'integer') {
      stats[ rownames(rv), colnames(rv)] = stats[rownames(rv), colnames(rv)] + rv
    }
  }
    
  stats
}

.sch_fends_pair_orientation <- function(nm, fends, mindist=100, maxdist=1e4, chroms=sch_chroms, log_step=0.125)
{
  x = gextract(nm, intervals=gintervals.2d(chroms), band=c(-maxdist, -mindist), colnames='count')

  fends = fends[is.element(fends$chrom, paste0("chr", chroms)), ]
  fends$dir = (1:nrow(fends)) %% 2
  
  m = merge(x, fends, by.x=c('chrom1', 'start1'), by.y=c('chrom', 'start'))
  m = merge(m, fends, by.x=c('chrom2', 'start2'), by.y=c('chrom', 'start'))
  m$dist_b = paste0('b', log_step * floor(log2(m$start2 - m$start1)/log_step))
  table(paste0('f', m$dir.x, m$dir.y), m$dist_b)
}


################
#
count_unique_fends_per_cell <- function()
{
  nms = gtrack.ls(sch_track_base)
  
  commands = paste0(".nm_count_unique_fends_per_cell(nm=\"", nms, "\")", collapse=",")
  res = eval(parse(text=paste("gcluster.run(",commands,", max.jobs=300)")))
  stopifnot(any(lapply(res, function(r)  typeof(r$retv) == "integer"  )))

  sch_fends_per_cell <<- do.call("c", lapply(res, function(r) r$retv))
  names(sch_fends_per_cell) = nms

  write.table(as.data.frame(sch_fends_per_cell), sprintf("%s/sch_n_distinct_fends_per_cell.txt", sch_table_dir), quote=F, sep="\t")
  
}

.nm_count_unique_fends_per_cell <- function(nm)
{
  library(dplyr)
  
  x = gextract(nm, intervals=gintervals.2d.all())
  x = filter(x, as.character(chrom1) < as.character(chrom2) | as.character(chrom1) == as.character(chrom2) & start1 < start2)
  length(unique(paste(c(x$chrom1, x$chrom2), c(x$start1, x$start2))))
}

#===============================================================================
# New way to count AA/AB/BB cis counts per cell
#===============================================================================

calc_ab_contacts <- function(mindist=2e6, chroms=sch_chroms)
{
  message("calc_ab_contacts")
  d = .sch_get_pool_domains()
  ab = read.table(paste0(sch_ab_cluster_fn_pref, ".tab.kclust"), header=T, stringsAsFactors=T)
  rownames(ab) = ab$id
  ab = ab[d$id, ]
  d$ab = ab$clust

  nms = gtrack.ls(sch_track_base)

  commands = paste0(".calc_cell_ab_contacts(nm=\"", nms, "\", d=d, mindist=mindist, chroms=chroms)", collapse=",")
  res = eval(parse(text=paste("gcluster.run(",commands,")")))
  stopifnot(sum(unlist(lapply(res, function(r) { typeof(r$retv) == 'list'} ))) == length(res))
  counts = NULL

  for (i in 1:length(nms)) {
    message(i)
    counts = rbind(counts, res[[i]]$retv)
  }
  rownames(counts) = nms
  write.table(counts, sprintf("%s/cell_ab_counts.txt", sch_table_dir), quote=F, sep="\t")
  sch_ab_counts <<- counts
}

.calc_cell_ab_contacts <- function(nm, d, mindist=2e6, chroms=sch_chroms)
{
  v_cis = v_trans = data.frame(d1=c(0, 0, 1), d2=c(0, 1, 1))

  x_cis = gextract(nm, gintervals.2d(chroms), colnames='cov', band=c(-max(gintervals.all()$end), -mindist))

  tr_i = gintervals.2d.all()
  c1 = as.character(tr_i$chrom1)
  c2 = as.character(tr_i$chrom2)
  tr_i = tr_i[c1 < c2 & is.element(c1, paste0('chr', chroms)) & is.element(c2, paste0('chr', chroms)),]

  x_trans = gextract(nm, intervals=tr_i, colnames='cov')

  x = rbind(x_cis, x_trans)

  if (!is.null(x)) {
    x1 = data.frame(chrom=x$chrom1, start=x$start1, end=x$end1, id=1:nrow(x))
    x2 = data.frame(chrom=x$chrom2, start=x$start2, end=x$end2, id=1:nrow(x))
    t1 = gintervals.neighbors(x1, d, maxneighbors=1, mindist=0, maxdist=0)
    t2 = gintervals.neighbors(x2, d, maxneighbors=1, mindist=0, maxdist=0)

    ids = as.character(intersect(t1$id, t2$id))

    if (length(ids) > 0) {
      rownames(t1) = t1$id
      rownames(t2) = t2$id

      colnames(t1) = c(paste0(colnames(x1), 1), paste0(colnames(d), 2), 'dist')
      colnames(t2) = c(paste0(colnames(x2), 1), paste0(colnames(d), 2), 'dist')

      t12 = data.frame(cis=as.character(t1[as.character(ids), 'chrom1']) == as.character(t2[as.character(ids), 'chrom1']), dist=t2[as.character(ids), 'start2'] - t1[as.character(ids), 'end2'], d1=pmin(t1[as.character(ids), 'ab2'], t2[as.character(ids), 'ab2']), d2=pmax(t1[as.character(ids), 'ab2'], t2[as.character(ids), 'ab2']))

      if (sum(t12$cis & t12$dist > mindist) > 1) {
        v_cis = rbind(v_cis, t12[t12$cis & t12$dist > mindist, c('d1', 'd2')])
      }

      if (sum(!t12$cis) > 1) {
        v_trans = rbind(v_trans, t12[!t12$cis, c('d1', 'd2')])
      }
    }
  }

  cis_tab = table(v_cis) - 1
  trans_tab = table(v_trans) - 1

  data.frame(cis_AA=cis_tab["1", "1"], cis_AB=cis_tab["0", "1"], cis_BB=cis_tab["0", "0"], trans_AA=trans_tab["1", "1"], trans_AB=trans_tab["0", "1"], trans_BB=trans_tab["0", "0"], total=ifelse(is.null(x), 0, nrow(x)))
}
