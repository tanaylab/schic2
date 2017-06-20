source("analyzeScHiC.r")
source("plots.r")
source("data_generation.r")
source("compartments.r")
library("plotrix")
library("ggplot2")
library("shaman")

options(gmax.data.size=1e8)
options(gmultitasking=F)


# params for paper
plot_type <<- "ps"
rel_odir <<- "paper"
fig_factor <<- 2.3
#fres <<- 120
fres <<- 300
pointsize <<- 7

# params for poster
#rel_odir <<- "poster"
#fig_factor <<- 2
#fres <<- NA
#pointsize <<- 12




#system(sprintf("mkdir -p %s/%s", sch_fig_dir, rel_odir))
build_stats <- function()
{
  import_rna_seq_fastq()
  
  create_pooled_tracks()
  sch_gen_group_pools()
  sch_mit_align_test(nms=rownames(sch_batch))
  km = fig_s_domains_trans_cluster()
  calc_ab_contacts()
  wrapper_sch_calc_tad_exact_borders_insu_ds_trian_discard_unear()
  sch_create_domain_marg_cov()
  wrapper_sch_calc_dixon_borders_insu_ds_trian()

  # pooled shuffled and knn-score tracks
  options(shaman.sge_support=1)
  shaman_shuffle_hic_track(track_db=sch_groot, obs_track_nm=pool_tn, work_dir=sprintf("%s/", sch_rdata_dir))
  shaman_score_hic_track(track_db=sch_groot,work_dir=sprintf("%s/", sch_rdata_dir), obs_track_nms=pool_tn, score_track_nm=paste0(pool_tn, "_score"))
  shaman_score_hic_track(track_db=sch_groot,work_dir=sprintf("%s/", sch_rdata_dir), obs_track_nms=pool_tn, score_track_nm=paste0(pool_tn, "_lowerres_score"), near_cis=1e7)

  gen_post_mitotic_insu_and_shuffle()
}

paper_figs <- function()
{
  fig_gw_cmap()

  fig1()
  fig2()
  fig3()

  exfig1()
  exfig2()
  exfig3()
  exfig4()
  
  exfig6()
  exfig7()
  exfig9_dip()

  dft  = compile_features_table()
  
  # load haploid settings
  sch_load_tables("~/schic_res/serum_2ic/serum_2i_271116_params.r")
  fig4()

  exfig8()
  exfig9_hap()

  
  # for ExFig 6d
  fig_s_phasing_by_batches()

  hft  = compile_features_table(haploid=T)
}

fig1 <- function(haploid=F)
{
  fig1_total_cov()
  fig1_trans_frac()
  fig1_phasing_metrics()
  fig1_decay_sorted_hm()
  fig1_cis_chrom_obs_per_group()
}

fig2 <- function()
{
  fig2_sorted_batches_on_a_wheel()
  fig2_metrics_by_cc()
  fig2_index_sort_metrics()
}

fig3 <- function()
{
  fig3_insu_ab_by_cc()
  fig3_insu_clusters()
  fig3_conv_ctcf_loops_by_cc()
  fig3_conv_ctcf_group_enr_plots()
  fig3_intra_domain_decay_by_tor()
}

# load the haploid settings before producing Fig4: 
#   sch_load_tables("~/schic_res/serum_2ic/serum_2i_271116_params.r")
fig4 <- function()
{
  fig1_decay_sorted_hm(haploid=T)
  r = fig_s_a_score_ideograms(tn=pool_tn, ref_tn="scell.nextera.pool_good_hyb_2i_idx_sort_es", lim=c(0.54, 0.78))
  r = plot_domains_weighted_a_score_across_groups()
  fig_s_pc_pc_contact_enrich(zlim=c(-1, 1), exp_zlim=c(-9.8, -7.8))
}


exfig1 <- function()
{
  fig_s_saturation()
  fig_s_qc()
  fig_s_fends_pairs_orientation()
}

exfig2 <- function()
{
  fig_s_trans_examples()
  fig_s_chrom_chrom_total_contacts()
}

exfig3 <- function()
{
  fig_s_pool_knn_norm_cis_maps()
  r = fig_s_a_score_ideograms()
  km = fig_s_domains_trans_cluster()
  fig_s_compare_tracks_insu()
}

exfig4 <- function()
{
  dec_refs = fig_s_group_decay_and_ensemble_ref()
  # ExtFig 4b cmaps are produced by fig_gw_cmap()
  fig_s_repli()
  fig_s_decay_and_repli_scatters()
  # NLDR projection is produced in python, see exfig_4e.py 
}

exfig5 <- function()
{
  fig_s_decay_kmeans()
}

exfig6 <- function()
{
  # 6a by fig_s_decay_kmeans()
  fig_s_2fold_cv_phasing()
  fig_s_chrom_metrics_phasing_by_cc()
  fig_s_phasing_by_batches()
  # 6g by fig2_index_sort_metrics()
  # 6h by fig_gw_cmap
}

exfig7 <- function()
{
  fig_s_compare_group_insu_trends()
  insu_s = fig_s_insu()
  fig_s_mitotic_ins()
  fig_s_trans_ab()
  fig_s_compare_conv_ctcf_loops_between_groups()
}

# load haploid params
exfig8 <- function()
{
  fig_s_saturation()
  fig_s_qc()

  fig2_metrics_by_cc(haploid=T)
  fig3_insu_ab_by_cc(haploid=T)
  
  s = fig3_conv_ctcf_loops_by_cc_haploid()
  fig3_intra_domain_decay_by_tor_haploid()
    
}

exfig9_hap <- function()
{
  # haploid
  fig_s_compare_tracks_a_score()
  r = plot_domains_weighted_a_score_across_groups(order_doms_by="kmean")

  fig_s_hap_dip_pc_contig_table()

  # 9e already produced in fig4()
  fig_s_domains_epi_stats()
}

exfig9_dip <- function()
{
  # diploid
  r = plot_domains_weighted_a_score_across_groups_diploid()
  fig_s_pc_pc_contact_enrich(zlim=c(-1, 1), exp_zlim=c(-9.8, -7.8))
  fig_s_hap_vs_dip_rna_seq_over_genes()
}

fig1_total_cov <- function(width=190, height=270, xlim=NULL, remove_nodig=T, good_cells=T)
{
  x = sch_chrom_stat

  if (good_cells) {
    x = x[sch_good_cells, ]
  }
  if (remove_nodig) {
    x = x[, grep("self", colnames(x), invert=T)]
  }
  v = rowSums(x)/2

  if (is.null(xlim)) {
    xlim = range(v)
  }
  fig_png(sprintf("%s/%s/Fig1_%s_%s_cov.png",sch_fig_dir, rel_odir, ifelse(good_cells, "good", "all"), ifelse(remove_nodig, "wo_nodig", "total")), width, height)
  plot(density(pmin(pmax(v, xlim[1]), xlim[2])), main="", xlab="coverage", ylab="density", lwd=2, xlim=xlim)
  print(quantile(v, c(0.25, 0.5, 0.75)))
  abline(v=median(v), col='blue', lty=2, lwd=2)
  #abline(v=sch_min_cov, col='red', lty=2, lwd=2)
  grid(col='darkgray')
  dev.off()
}

fig1_trans_frac <- function(width=190, height=270, good_cells=T)
{
  x = sch_glob_decay_res
  
  if (good_cells) {
    x = x[sch_good_cells,]
  }  
  v = x[, 'trans'] / rowSums(x)

  fig_png(sprintf("%s/%s/Fig1_trans_frac_%s.png",sch_fig_dir, rel_odir, ifelse(good_cells, "good", "all")), width, height)
  plot(density(pmin(v, sch_max_trans)), main="", xlab="%trans", ylab="density", lwd=2, xlim=c(0,sch_max_trans))
  print(median(v))
  abline(v=median(v), col='blue', lty=2, lwd=2)
  #abline(v=sch_max_trans, col='red', lty=2, lwd=2)
  grid(col='darkgray')
  dev.off()
}

fig1_decay_scatters <- function(width=180, height=220, xfield="high_near_f", yfield="far_tightness", xlim=NULL, ylim=c(0.03, 0.1), shift_to_lim=T, cex=0.3, col="black", col_by_cond=F, x=NULL, y=NULL, alpha=0.2, vlines=NULL, vcol='red', vlty=2, vlwd=2, data=sch_decay_metrics)
{

  if (is.null(x)) {
    x = data[, xfield]
  }

  if (is.null(y)) {
    y = data[, yfield]
  }

  if (is.null(xlim)) {
    xlim = range(x)
  }

  if (is.null(ylim)) {
    ylim = range(y)
  }

  if (shift_to_lim) {
    x[ x < xlim[1]] = xlim[1]
    x[ x > xlim[2]] = xlim[2]
    y[ y < ylim[1]] = ylim[1]
    y[ y > ylim[2]] = ylim[2]
  }

  if (col_by_cond) {
    col = rep("darkgray", length(x))
    col[data$cond %in% names(sch_cond_colors)] = unlist(sch_cond_colors[data$cond])
  }

  fig_png(sprintf("%s/%s/Fig1_dec_scat_%s_vs_%s%s.png",sch_fig_dir, rel_odir, xfield, yfield, ifelse(col_by_cond, "_by_cond", "")), width, height)
  plot(x, y, main="", xlab=xfield, ylab=yfield, pch=19, cex=cex, xlim=xlim, ylim=ylim, col=addalpha(col, alpha))
  if (col_by_cond) {
    ind = col != "darkgray"
    points(x[ind], y[ind], pch=19, cex=cex, col=addalpha(col[ind], alpha))
  }
  grid(col='gray', lwd=0.5)
  if (!is.null(vlines)) {
    abline(v=vlines, col=vcol, lty=vlty, lwd=vlwd)
  }

  dev.off()
}

fig1_phasing_metrics <- function (width=750, height=280, cex=0.7*fig_factor, alpha=1)
{
  m = sch_decay_metrics
  
  fig_png(sprintf("%s/%s/fig1_phasing_metrics_by_group.png", sch_fig_dir, rel_odir), width, height)
  par(mfrow=c(1,3))
  par(mar=c(3,3,1,1))
  
  plot(m$f_mitotic_band, m$f_near_band, pch=19, cex=cex, col=addalpha(mini_qualitative_colors[as.numeric(m$group)], alpha), xlab="%mitotic", ylab="%near", main="")
  grid(col='black', lwd=0.5)
  par(mar=c(3,0,1,4))
  plot(m$early_f/min(m$early_f), m$f_near_band,  pch=19, cex=cex, col=addalpha(mini_qualitative_colors[as.numeric(m$group)], alpha), yaxt='n', ylab="", xlab="repli-score", main="")
  grid(col='black', lwd=0.5)

  plot(m$far_mu, m$f_near_band,  pch=19, cex=cex, col=addalpha(mini_qualitative_colors[as.numeric(m$group)], alpha), yaxt='n', ylab="", xlab="far mu", main="")
  grid(col='black', lwd=0.5)

  dev.off()

}

fig1_decay_sorted_hm <- function(zlim=c(0, 0.025), colspec=seq_colspec, ncols=256, width=fig_factor*2600, height=fig_factor*800, add_ab=T, plot_lines=F, haploid=F, early_zlim=c(0.45, 0.75), ab_zlim=c(0.15, 0.95))
{
 # plot_legend(ofn=sprintf("%s/%s/Fig1_decay_hm_legend.png",sch_fig_dir, rel_odir), zlim=zlim, crp=colorRampPalette(colspec)(ncols), "fraction")

  vlines=NULL
  if (plot_lines) {
    vlines = cumsum(table(sch_decay_metrics$group)) / nrow(sch_decay_metrics)
  }
  
  op = plot_decay_mat_sorted_by_x(zlim=zlim, order_by="given", name="good_", cells_names_ord=rownames(sch_decay_metrics)[order(sch_decay_metrics$ord)], remove_trans=T, add_dist_grid=F, add_ab=add_ab, add_nodig_trans=F, remove_unear_dist=sch_remove_near_dist, width=width, height=round(ifelse(haploid, 1.33, 1) * height), ext_ofn=sprintf("%s/%s/Fig1_decay_hm_new%s_%s.png", sch_fig_dir, rel_odir, ifelse(plot_lines, "_glines", ""), ifelse(haploid, "hap", "dip")), vlines=vlines, ec_zlim=early_zlim, ab_zlim=ab_zlim, colspec=colspec)
}


fig1_cis_chrom_obs_per_group <- function(odir=sprintf("%s/%s/cmaps/", sch_fig_dir, rel_odir), chr="chr11", slice_size=60, bw=1e5, glob_dens_lim=c(1e-18, 5e-16), colspec=c(seq_colspec, 'black'))
{
  odir = sprintf("%s/%s_per_group", odir, chr)
  system(paste("mkdir -p", odir))
  
  mv = sch_decay_metrics

  for (i in sort(unique(as.numeric(mv$group)))) {
    cmv = mv[mv$group == i, ]
    from = pmax(1, round((nrow(cmv) - slice_size)/2))
    to   = pmin(nrow(cmv), round((nrow(cmv) + slice_size)/2))

    from = which(rownames(mv) == rownames(cmv)[from])
    to = which(rownames(mv) == rownames(cmv)[to])
    message(paste(i, from, to))

    .sch_plot_group_mean_mat(odir, rownames(mv)[from:to], scope=gintervals.2d(chr), bwx=bw, fn_pref=sprintf("gr%d_pool_%d-%d", i, from, to), type="kde", add_axis=T, rotate=F, force_size=-1, glob_dens_lim=glob_dens_lim, shades=colorRampPalette(colspec))  

  }
}

#############
fig2_sorted_batches_on_a_wheel <- function(width=800, height=800, conds=c("1CD_G1", "1CD_eS", "1CD_mS", "1CD_lS_G2"), names=c("G1", "Early S", "Mid S", "Late S/G2"), cols=c("blue", rgb(0.7, 0.3, 0.7) , rgb(0.48, 0.76, 0.48), rgb(0, 0.37, 0)), start_at=pi/2, line_len=0.7, name="Fig2")
{
  mv = sch_decay_metrics[ sch_decay_metrics$valid, ]
  mv = mv[order(mv$ord, decreasing=T), ]

  x = data.frame(cond=mv$cond, rad=start_at + (1:nrow(mv)) * 2 * pi / nrow(mv), layer=NA)

  x = x[is.element(x$cond, conds), ]
  x$layer = unlist(sapply(x$cond, function(y) { which (conds == y) }))

  r = 3 * length(conds) + 1

  fig_png(sprintf("%s/%s/%s_cond_wheel.png",sch_fig_dir, rel_odir, name), width, height)

  par(mar=c(1,1,1,1))
  plot.new()
  plot.window(xlim=c(-r, r), ylim=c(-r, r))
  y = 2 * length(conds)
  segments((y + x$layer ) * cos(x$rad), (y + x$layer) * sin(x$rad), (y + x$layer + line_len) * cos(x$rad), (y + x$layer + line_len) * sin(x$rad), col=cols[x$layer], lwd=0.5 * fig_factor)
  legend("center", legend=names, lty=1, col=cols, bty='n', cex=1.5*fig_factor, lwd=2*fig_factor, seg.len=line_len)
  draw.circle(c(0, 0), c(0, 0), c(y + line_len, r), nv=360, border="darkgray")
  segments((y - 2) * cos(start_at), (y - 2) * sin(start_at), y * cos(start_at), y * sin(start_at), lwd=0.5 * fig_factor)
  arrows((y - 1) * cos(start_at), (y - 1) * sin(start_at), (y - 1) * cos(start_at) + 1, (y - 1) * sin(start_at), code=2, length=0.05*fig_factor)
  dev.off()

}

fig2_sorted_batches_on_a_wheel_haploid <- function()
{
  fig2_sorted_batches_on_a_wheel(conds=names(sch_cond_colors), names=c("None", "G1", "G1/S"), cols=unlist(sch_cond_colors), start_at=pi/2, line_len=0.7)
}

fig_s_all_sorted_batches_on_a_wheel <- function()
{
  nms = unique(names(sch_cond_colors))
  cols = unlist(unique(sch_cond_colors))
  fig2_sorted_batches_on_a_wheel(conds=nms, names=nms, cols=cols, start_at=pi/2, line_len=0.7, name="FigS_all")
}

fig_s_serum_sorted_batches_on_a_wheel <- function()
{
  nms = c('1CDS1', '1CDS2')
  cols = c("dodgerblue1", "orange")
  fig2_sorted_batches_on_a_wheel(conds=nms, names=nms, cols=cols, start_at=pi/2, line_len=0.7, name="FigS_serum")
}

fig_s_phasing_by_batches <- function()
{
  m = sch_decay_metrics

  m$batch_id = sch_batch[rownames(m), 'batch_id']
  
  b = as.data.frame(unique(cbind(as.numeric(m$batch), m$batch_id, m$cond)))
  colnames(b) = c('batch', 'batch_id', 'cond')
  b = b[order(as.numeric(as.character(b$batch))), ]
  b$n = 1:nrow(b)
  rownames(b) = as.character(b$batch_id)

  m$batch_n = b[as.character(m$batch_id), 'n']
  
  fig_png(sprintf("%s/%s/fig_s_phasing_by_batches.png", sch_fig_dir, rel_odir), 720, 300)
  plot.new()
  par(mar=c(3,5,1,5))

  plot.window(xlim=range(m$ord), ylim=c(0, nrow(b)))

  rect(1, 0, max(m$ord), nrow(b), border='black', col='white', lwd=2)
  segments(m$ord, m$batch_n - 1, m$ord, m$batch_n, col=unlist(sch_cond_colors[m$cond]))
  axis(1)
  axis(2, at=seq(1, nrow(b)) - 0.5, labels=rownames(b), las=2, tick=F, hadj=1, line=-2)
  axis(4, at=seq(1, nrow(b)) - 0.5, labels=b$cond, las=2, tick=F, hadj=0, line=-2)

  segments(cumsum(table(m$group)), 1, cumsum(table(m$group)), nrow(b))
  
  dev.off()
  
}

##
fig2_metrics_by_cc <- function(width=fig_factor*360, height_per_panel=fig_factor*110, sw=101, cex=fig_factor*0.4, lwd=fig_factor*2, collapse_q=c(0.005, 0.995), col_by_cond=T, only_facs=F, haploid=F, cols=NULL)
{
  mv = sch_decay_metrics[ sch_decay_metrics$valid, ]
  
  align = read.table(sprintf("%s/trans_contacts_corr.tab", sch_table_dir), header=T, stringsAsFactors=F)
  rownames(align) = align$nm
  align = align[rownames(mv),]

  #cov = sapply(rownames(mv), function(nm) { message(nm); gsummary(nm)["Sum"]/2 })
  cov = rowSums(sch_chrom_stat[rownames(mv),])/2
  #reads = unique(read.table(sprintf("%s/reads_per_cell.txt", sch_table_dir)))
  #colnames(reads) = c('n', 'nm')
  #rownames(reads) = paste0('scell.nextera.', as.character(reads$nm))
  #reads = reads[rownames(mv),]
  x = data.frame(early_cov=mv$early_f/min(mv$early_f), cov=cov, trans_f=align$ntrans/align$n,  tr_align=align$c_align)

  if (is.null(cols)) {
    cols = get_cond_cols(haploid)
  }

  if (only_facs) {
    x = x[grep("1CD_", mv$cond), ]
    cols = cols[grep("1CD_", mv$cond)]
    
  }
  .plot_points_and_trend(x, ofn=sprintf("%s/%s/Fig2_metrics_by_cc_sw%d.png",sch_fig_dir, rel_odir, sw), cols, width=width, height_per_panel=height_per_panel, sw=sw, col_by_cond=col_by_cond, collapse_q=collapse_q, lcol='black', cex=cex, lwd=lwd, outlier_col=ifelse(haploid, 'black', 'red'))
}

##
fig_gw_cmap <- function(chroms=sch_chroms, binsize=8e6, zlim=NULL, colspec = c("white", "blue", "darkblue", "red",  "black", "yellow", "purple"))
{
  odir = sprintf("%s/%s/cmaps/gw/bs_%s_%g_%g",sch_fig_dir, rel_odir, n2str(binsize, 0), ifelse(is.null(zlim), NA, zlim[1]), ifelse(is.null(zlim), NA, zlim[2]))
  mv = sch_decay_metrics[sch_decay_metrics$valid, ]

  system(paste("mkdir -p", odir))

  ofns = sapply(1:nrow(mv), function (i) { sprintf("%s/%04d_g%s_%s.png", odir, i, mv[i, 'group'], gsub("scell.nectera.", "", rownames(mv)[i])) })
  
  commands = paste0("sch_plot_genome_wide(\"", rownames(mv), "\", \"", ofns, "\", cis_add_cols=colspec, trans_cols=c(rep('white', 3), 'blue'), chroms=chroms, binsize=binsize, zlim=zlim)", collapse=",")
  print(commands)
  res = eval(parse(text=paste("gcluster.run(",commands,")")))
  res
}

####
fig3_insu_ab_by_cc <- function(width=fig_factor*220, height_per_panel=fig_factor*150, sw=101, cex=fig_factor*0.4, lwd=fig_factor*2, collapse_q=c(0.005, 0.995), col_by_cond=T, insu_data=NULL, rebuild=F, filter_cov_domains=F, plot_trans_ab=F, haploid=F)
{
  mv = sch_decay_metrics[ sch_decay_metrics$valid, ]

  cis_ab = calc_ab_enrich(sch_ab_counts[rownames(mv),  grep("cis", colnames(sch_ab_counts))], "cis")
  trans_ab = calc_ab_enrich(sch_ab_counts[rownames(mv),  grep("trans", colnames(sch_ab_counts))], "trans")

  if (is.null(insu_data)) {
    insu_name = sprintf("_tad_borders_gt%s%s", n2str(sch_remove_near_dist), ifelse(filter_cov_domains, "_good_cov", ""))
    insu_data = insu_by_cc_grouped(sch_decay_metrics, insu_name=insu_name)
  }

  x = data.frame(insu=-insu_data$mean_ins_per_cell,
    cis_ab_comp=cis_ab$deplAB,
    trans_ab_comp=trans_ab$deplAB)

  cols = get_cond_cols(haploid)
  
  ofn_base = sprintf("%s/%s/Fig3_insu_ab_by_cc_sw%d%s",sch_fig_dir, rel_odir, sw, ifelse(filter_cov_domains, "_covD", ""))
  if (plot_trans_ab) {
    xa = x[,1:3]
  }
  else {
    xa = x[,1:2]
  }
  .plot_points_and_trend(xa, ofn=paste0(ofn_base, "_a.png"), cols, width=width, height_per_panel=height_per_panel, sw=sw, col_by_cond=col_by_cond, collapse_q=collapse_q, lcol='black', vertical=F, cex=cex, lwd=lwd, show_yaxis=T, outlier_col=ifelse(haploid, 'black', 'red'))


  insu_data
}

# fig3_conv_ctcf_loops_by_cc submitted: function(width=fig_factor*220, height_per_panel=fig_factor*150, sw=101, cex=0.4, collapse_q=c(0.005, 0.995), col_by_cond=T, rebuild=F, ctcf_chip_q=0.99, tor_horiz=2e5, band=c(2e5, 1e6), conv_ctcf=NULL, n_contact_scale=1e6, loop_widths=c(1e4, 3e4), conv_ctcf_pool_q=0, group_size_for_boxplot=150, boxplot_comp_inds=matrix(c(1,4,4,10,10,14, 14, 17), 4, 2, byrow=T), show_pvals_as_nums=T, min_d_score=60, score_tn=paste0(pool_tn, "_score"), en_boxplot_ylim=c(0, 1.65), haploid=F)
fig3_conv_ctcf_loops_by_cc <- function(width=fig_factor*220, height_per_panel=fig_factor*150, sw=101, lwd=fig_factor*2, cex=fig_factor*0.4, collapse_q=c(0.005, 0.995), col_by_cond=T, rebuild=F, ctcf_chip_q=0.99, tor_horiz=2e5, band=c(2e5, 1e6), conv_ctcf=NULL, n_contact_scale=1e6, loop_widths=c(1e4, 3e4), conv_ctcf_pool_q=0, group_size_for_boxplot=100, boxplot_comp_inds=matrix(c(1,4,4,7,7,10,10,13, 13, 17), 5, 2, byrow=T), show_pvals_as_nums=T, min_d_score=60, score_tn=paste0(pool_tn, "_score"), en_boxplot_ylim=c(0, 1.8), haploid=F, loop_foc_enr_q=0.8, count_contact_once=T, external_loops=NULL, fn_pref=sprintf("%s/%s/fig3_loops_",sch_fig_dir, rel_odir), name="")
{
  plot_conv_ctcf_counts(width=width, height_per_panel=height_per_panel, sw=sw, cex=cex, collapse_q=collapse_q, col_by_cond=col_by_cond, rebuild=rebuild, ctcf_chip_q=ctcf_chip_q, tor_horiz=tor_horiz, band=band, conv_ctcf=conv_ctcf, n_contact_scale=n_contact_scale, loop_widths=loop_widths, loop_foc_enr_q=loop_foc_enr_q, conv_ctcf_pool_q=conv_ctcf_pool_q, group_size_for_boxplot=group_size_for_boxplot, boxplot_comp_inds=boxplot_comp_inds, show_pvals_as_nums=show_pvals_as_nums, fn_pref=fn_pref, plot_horiz=T, min_d_score=min_d_score, score_tn=score_tn, en_boxplot_ylim=en_boxplot_ylim, haploid=haploid, count_contact_once=count_contact_once, external_loops=external_loops, name=name, lwd=lwd)
}

fig3_conv_ctcf_loops_by_cc_haploid <- function(rebuild=F, conv_ctcf=NULL)
{
  fig3_conv_ctcf_loops_by_cc(rebuild=rebuild, conv_ctcf=conv_ctcf, group_size_for_boxplot=70, boxplot_comp_inds=matrix(c(1,7,7,13,13,18), 3, 2, byrow=T), show_pvals_as_nums=T, score_tn=paste0(pool_tn, "_score"), collapse_q=c(0.001, 0.999), haploid=T, en_boxplot_ylim=c(-0.5, 2.8))
}

fig_s_rao_loops_enr_by_cc <- function(band=c(1e5, 1e6), rebuild=F)
{
  rao = read.table("~/schic_res/hyb_mm9_2i/tables/rao2014/GSE63525_CH12-LX_HiCCUPS_looplist.txt", header=T, sep="\t")
  rao_l = gintervals.2d(rao$chr1, rao$x1, rao$x2, rao$chr2, rao$y1, rao$y2)

  fig3_conv_ctcf_loops_by_cc(external_loops=rao_l, fn_pref=sprintf("%s/%s/fig_s_rao_loops_",sch_fig_dir, rel_odir), band=band, name="rao_", rebuild=rebuild)
}

fig3_conv_ctcf_group_enr_plots <- function (conv_ctcf=NULL, ctcf_chip_q=0.99, min_score=60, extend_by=50e3, nbins=c(7, 31), zlim=c(-1, 1), colspec=diver_colspec)

{
  #colspec = colspec=c('darkblue', 'blue', 'white', 'red', 'darkred')
  orig_multi = options("gmultitasking")
  options(gmultitasking=T)

  for (bins in nbins) {
    conv_ctcf = plot_contacts_between_score_filtered_conv_ctcf_by_pooled_group(extend_by=extend_by, nbins=bins, conv_ctcf=conv_ctcf, ctcf_chip_q=ctcf_chip_q, min_score=min_score, odir=sprintf("%s/%s",sch_fig_dir, rel_odir), colspec=colspec, zlim=zlim)
  }
  options(gmultitasking=orig_multi)

  conv_ctcf
}

fig3_insu_clusters <- function(insu_data=NULL, width=330, height_per_panel=165, sw=101, collapse_q=c(0.005, 0.995), cex=fig_factor*0.2, lwd=fig_factor*2, pcol=rev(seq_colspec)[3], colspec=diver_colspec, single_tor_col=T, group_size=100)
{

  if (is.null(insu_data)) {
    insu_data = insu_by_cc_grouped(sch_decay_metrics, group_size=group_size)
  }

  km = insu_data$km
  bt = insu_data$bord_tor
  bt = bt[names(km$cluster), ]

  cord = order(tapply(bt$tor, bt$cluster, mean, na.rm=T))

  x = -t(insu_data$yc)
  x = x[, cord]

  ic = split(km$cluster, km$cluster)

  k = length(km$size)
  fig_png(sprintf("%s/%s/Fig3_insu_clus_tor_scatters.png",sch_fig_dir, rel_odir), width=height_per_panel*1.2, height=height_per_panel * k)
  par(mar=c(1,4,0,1))
  layout(matrix(1:k, k, 1))


  #xlim = range(bt$tor, na.rm=T)
  #ylim = range(abs(bt$up_tor - bt$down_tor), na.rm=T)
  xlim = range(bt$down_tor, na.rm=T)
  ylim = range(bt$up_tor, na.rm=T)

  for (i in seq_along(cord)) {
    cbt = bt[ names(ic[[cord[i]]]),]

    message(paste(i, nrow(cbt), mean(cbt$tor, na.rm=T)))
    #plot(cbt$tor, abs(cbt$up_tor - cbt$down_tor), pch=19, cex=0.2, xlab="", ylab=paste0("c", i), xaxt='n', xlim=xlim, ylim=ylim)
    if (single_tor_col) {
      col = darkblue_col
    }
    else {
      col = colorRampPalette(colspec)(256)[vals_to_cols(cbt$tor, c(-1.5,0, 1.5), 128)]
    }
    plot(cbt$up_tor, cbt$down_tor, pch=19, cex=cex, xlab="", ylab=paste0("c", i), xaxt='n', xlim=xlim, ylim=ylim, asp=1, col=col)
    grid(col='darkgray', lwd=0.5*fig_factor)
  }
  axis(1, cex=0.5)
  dev.off()

  colnames(x) = paste0("c", 1:ncol(x))

  print(apply(x, 2, which.max))
  .plot_points_and_trend(x, ofn=sprintf("%s/%s/Fig3_insu_clus_by_cc_sw%d.png",sch_fig_dir, rel_odir, sw), cols=NULL, pcol=pcol, width=width, height_per_panel=height_per_panel, sw=sw, col_by_cond=F, collapse_q=collapse_q, glob_trend=cyclic_rollmean(-insu_data$mean_ins_per_cell, sw), glob_col="black", glob_ylim=quantile(x, collapse_q), lcol=rev(seq_colspec)[1], cex=cex, lwd=lwd)


  insu_data
}

fig3_intra_domain_decay_by_tor <- function(dom_dist_cutoffs=c(0, 17.9, 18.9, 19.9, 20.5), decay_step=0.1, chroms=sch_chroms, dom_cutoff_q=0.15, width_per_panel=fig_factor*180, height_per_panel=fig_factor*140, sw=101, cex=fig_factor*0.4, lwd=fig_factor*2, collapse_q=c(0.002, 0.998), col_by_cond=T, ylim_qs = c(0.005, 0.995), rebuild=F, group_cols=c('#74a9cf', '#023858', 'red', 'green', 'darkgreen'), plot_decay_by_group=F, n_slices=5, xlims_ends=c(17.5, 18, 18.5, 19), ylims_starts=c(0.015, 0.0175, 0.0175, 0.015), group_size_for_boxplot=100, haploid=F)
{

  intra = gen_domains_decay_by_size_and_tor(dom_dist_cutoffs=dom_dist_cutoffs, decay_step=decay_step, chroms=chroms, dom_cutoff_q=dom_cutoff_q, rebuild=rebuild)

  tor_types = c("EE", "E", "L", "LL")

  mv = sch_decay_metrics[sch_decay_metrics$valid, ]

  cols = get_cond_cols(haploid)

  for (i in 1:(length(dom_dist_cutoffs)-1)) {
    df = NULL
    cylim = NULL

    y = list()
    for (tor_s in tor_types) {
      x = intra[[i]][[tor_s]]
      x_n = x / rowSums(x)
      x_n[rowSums(x) == 0, ] = 0
      v = as.matrix(x_n) %*% 2^as.numeric(colnames(x_n))
      df = rbind(df, as.vector(v))


      if (plot_decay_by_group) {
        y[[tor_s]] = do.call("rbind", lapply(split(1:nrow(x_n), mv$group), function(ind) { colMeans(x_n[ind, ]) }))
      }
      else {
        y[[tor_s]] = do.call("rbind", lapply(split(1:nrow(x_n), ceiling(n_slices * ((1:nrow(x_n))/nrow(x_n)))), function(ind) { colMeans(x_n[ind, ]) }))
      }
      cylim = range(c(cylim, range(y[[tor_s]])))
    }

    fig_png(sprintf("%s/%s/FigX_intra_domain_decay_by_%s_d%d_%s-%s_sw%d.png",sch_fig_dir, rel_odir, ifelse(plot_decay_by_group, "group", sprintf("%d_slices", n_slices)), i, n2str(2**dom_dist_cutoffs[i], 1), n2str(2**dom_dist_cutoffs[i+1], 1), sw), 120 * length(tor_types), 120)
    layout(matrix(1:length(y), 1, length(y)))
    par(mar=c(2,2,1,1))
    for (tor_s in tor_types) {
      cy = y[[tor_s]]
      ylim = c(ylims_starts[i], cylim[2])
      xlim = c(min(as.numeric(colnames(cy))), xlims_ends[i])
      plot(as.numeric(colnames(cy)), cy[1,], col=NA, main=paste(i, tor_s), xlab="distance (log2)", ylab="fraction", type='l', xlim=xlim, ylim=ylim)
      apply(cbind(group_cols, cy), 1, function(x) { lines(as.numeric(colnames(cy)), x[-1], col=x[1], lwd=1) } )
      grid(col='darkgray')
      if (plot_decay_by_group) {
        legend("bottomleft", legend=paste0('g', rownames(cy)), col=group_cols, lwd=1, bty='n', cex=0.8, ncol=2)
      }
    }
    dev.off()

    rownames(df) = tor_types

    message(sprintf("%d\t%s\t%s", i, n2str(quantile(df, ylim_qs[1]), 2), n2str(quantile(df, ylim_qs[2]), 2)))

    .plot_points_and_trend(t(df), ofn=sprintf("%s/%s/FigX_intra_domain_decay_by_cc_d%d_%s-%s_sw%d.png",sch_fig_dir, rel_odir, i, n2str(2**dom_dist_cutoffs[i], 1), n2str(2**dom_dist_cutoffs[i+1], 1), sw), cols, width_per_panel=width_per_panel, height_per_panel=height_per_panel, sw=sw, col_by_cond=col_by_cond, collapse_q=collapse_q, lcol='black', glob_ylim=quantile(df, ylim_qs), vertical=F, group_size_for_boxplot=group_size_for_boxplot, cex=cex, lwd=lwd, outlier_col=ifelse(haploid, 'black', 'red'))
  }

  intra
}

fig3_intra_domain_decay_by_tor_haploid <- function(rebuild=F)
{
  fig3_intra_domain_decay_by_tor(cond_cols=sch_cond_colors, rebuild=rebuild, width_per_panel=fig_factor*205, height_per_panel=fig_factor*160)
}

#fig_s_compare_group_insu_trends = function(  g_tns=c("scell.nextera.pool_group_2_g1_diploid", "scell.nextera.pool_group_3_early_s_diploid", "scell.nextera.pool_group_4_mid_s_g2_diploid"), g_nms=c("G1", "early-S", "mid-S/G2"), chroms=c(1:5), xlim=c(32e6, 38e6), ylim_q=c(0.001, 0.999), gcols=c('blue', 'red', 'darkgreen'), height_per_chrom=60, width=280, scope=NULL)
fig_s_compare_group_insu_trends <- function(  g_nms = c('post_m', 'g1', 'early_s', 'mid_s_g2', 'pre_m'), chroms=c(1:5), xlim=c(32e6, 38e6), ylim_q=c(0.001, 0.999), gcols=mini_qualitative_colors[1:5], height_per_chrom=60, width=280, scope=NULL, res=1000, selected_groups=2:4)
{
  g_nms = g_nms[selected_groups]
  g_tns=paste(pool_tn, "group", selected_groups, g_nms, sep="_")

  g_ins = paste0(g_tns, "_ins_", format(ins_scale, scientific=F), "s")
  if (!gtrack.exists(g_ins[1])) {
    for (i in seq_along(g_ins)) {
      gtrack.2d.gen_insu_track(g_tns[[i]], ins_scale, res=res, new_track=g_ins[[i]])
    }
  }
  
  insp = gextract(g_ins, intervals=gintervals(chroms, xlim[1], xlim[2]), iterator=res, colnames=g_nms)
  fig_png(sprintf("%s/%s/S_ins_trend_by_grp_examples.png",sch_fig_dir, rel_odir), width=width, height=height_per_chrom*length(chroms))
  par(mfrow=c(length(chroms), 1))
  par(mar=c(0.5,4,0.5,0.5))

  #ylim = range
  ylim = quantile(-insp[, g_nms], ylim_q, na.rm=T)

  for (i in seq_along(chroms)) {
    cg = insp[insp$intervalID == i, ]
    plot(cg$start, -cg[, g_nms[1]], type='l', col=NA, ylim=ylim, ylab=sub("chr", "", chroms[i]), main="", xaxt='n')
    grid(col='black', lwd=0.5)
    if (i == 1) {
      legend("topright", legend=g_nms, col=gcols, lty=1, bty='n', cex=0.8, ncol=length(g_nms))
    }
    apply(rbind(gcols, -cg[, g_nms]), 2, function(x) { lines(cg$start, x[-1], col=x[1]) })
  }
  dev.off()

  gb = get_tad_borders_per_group(g_nms=g_nms, g_tns=g_tns, scope=scope)
  gvn = -gb$b_ins_vals[, seq(4, by=1, length=length(g_ins))]
  g_ins_mean = lapply(g_ins, function(s) { -gsummary(s)[['Mean']] })
  names(g_ins_mean) = colnames(gvn)
  o = combn(length(g_ins), 2)

  fig_png(sprintf("%s/%s/S_ins_by_grp_scatters.png",sch_fig_dir, rel_odir), width=400, height=150)
  par(mfrow=c(1, length(g_ins)))
  par(mar=c(4,4,0.5,0.5))

  lim = c(min(unlist(g_ins_mean)), quantile(gvn, 0.99, na.rm=T))

  for (i in 1:ncol(o)) {
    xs = colnames(gvn)[o[1,i]]
    ys = colnames(gvn)[o[2,i]]
    message(paste(xs, ys))
    plot(gvn[,xs], gvn[,ys], xlab=xs, ylab=ys, main="", col=addalpha(darkblue_col, 0.2), pch=19, cex=0.3, xlim=lim, ylim=lim, asp=1)
    abline(h=g_ins_mean[[ys]], col='black', lty=2)
    abline(v=g_ins_mean[[xs]], col='black', lty=2)
    abline(a=0, b=1, col='black')
  }
  dev.off()
}

#fig_s_compare_conv_ctcf_loops_between_groups = function(g_tns=c("scell.nextera.pool_group_2_g1_diploid", "scell.nextera.pool_group_3_early_s_diploid", "scell.nextera.pool_group_4_mid_s_g2_diploid"), g_nms=c("G1", "early-S", "mid-S/G2"), chroms=sch_chroms, gcols=c('blue', 'red', 'darkgreen'), size=220, min_d_score=60, score_tn="scell.nextera.pool_good_hyb_es_b_score_k100", conv_ctcf=NULL, enr=NULL, pool_loop_extends=c(1e4, 3e4), ctcf_chip_q=0.99, ex_chroms=c(1:5), xlim=c(32e6, 38e6))
fig_s_compare_conv_ctcf_loops_between_groups <- function(g_nms = c('post_m', 'g1', 'early_s', 'mid_s_g2', 'pre_m'), chroms=sch_chroms, gcols=mini_qualitative_colors[1:5], size=220, min_d_score=60, score_tn=paste0(pool_tn, "_score"), conv_ctcf=NULL, enr=NULL, point_size=0.05, loop_size=4, loop_stroke=1, pool_loop_extends=c(1e4, 3e4), ctcf_chip_q=0.99, ex_chroms=c(1:5), xlim=c(26e6, 44e6), selected_groups=2:4, plot_knn_maps=T, plot_domains=F, alpha=0.5, cex=0.5*fig_factor, min_a_contacts=2)
{
  if (is.null(conv_ctcf)) {
    conv_ctcf = .get_conv_ctcf_with_tor(pool_q_cutoff=0, pool_loop_extend=pool_loop_extends[1], min_d_score=min_d_score, score_tn=score_tn, ctcf_chip_q=ctcf_chip_q)
  }

  if (plot_knn_maps) {
    for (chr in ex_chroms) {
      message(paste("plotting ", chr))
      knn_data = plot_knn_norm_cis_map(chr=chr, coords=xlim, bp_per_px=10e3, conv_ctcf=conv_ctcf, plot_loops=T, loop_size=loop_size, point_size=point_size, loop_stroke=loop_stroke, plot_also_domains=plot_domains, band=c(-2e6, 0))
      
    }
  }
  
  g_nms = g_nms[selected_groups]
  g_tns=paste(pool_tn, "group", selected_groups, g_nms, sep="_")
  
  if (is.null(enr)) {

    for (i in seq_along(g_tns))  {
      message(paste("extracting loop counts for", g_nms[i]))
      gvtrack.create("tn1", g_tns[i], "weighted.sum")
      gvtrack.create("tn2", g_tns[i], "weighted.sum")
      gvtrack.iterator.2d("tn1", -pool_loop_extends[1], pool_loop_extends[1], -pool_loop_extends[1], pool_loop_extends[1])
      gvtrack.iterator.2d("tn2", -pool_loop_extends[2], pool_loop_extends[2], -pool_loop_extends[2], pool_loop_extends[2])


      lc = gextract("ifelse(is.na(tn1), 0, tn1)", "ifelse(is.na(tn2), 0, tn2)", intervals=conv_ctcf, iterator=conv_ctcf, colnames=c('a', 'b'))
      lc[lc$a < min_a_contacts, 'a'] = NA
      enr = cbind(enr, lc$a / lc$b)
    }
    colnames(enr) = g_nms
  }

  o = combn(length(g_tns), 2)

  fig_png(sprintf("%s/%s/S_conv_ctcf_min%d_by_grp_scatters.png",sch_fig_dir, rel_odir, min_d_score), width=400,  height=150)
  par(mfrow=c(1, length(g_tns)))
  par(mar=c(4,4,0.5,0.5))

  
  unif_enr = (pool_loop_extends[1] / pool_loop_extends[2])**2
  enr_lr = log2(enr / unif_enr)
  lim = range(enr_lr, na.rm=T)
  #lim = c(min(enr_lr, na.rm=T), quantile(enr_lr, 0.99, na.rm=T))
  for (i in 1:ncol(o)) {
    xs = colnames(enr)[o[1,i]]
    ys = colnames(enr)[o[2,i]]
    message(paste(xs, ys))
    plot(enr_lr[,xs], enr_lr[,ys], xlab=xs, ylab=ys, main="", col=addalpha(darkblue_col, alpha), pch=19, cex=cex, xlim=lim, ylim=lim, asp=1, cex.lab=1.3)
    grid(col='black', lwd=0.5*fig_factor)
    #abline(a=0, b=1, lwd=1, col='red')
    abline(h=0, col='black', lwd=1*fig_factor)
    abline(v=0, col='black', lwd=1*fig_factor)
    legend("topright", legend=sprintf("r=%.2f", cor(enr_lr[,xs], enr_lr[,ys], use="comp")), col=NA, bty='n', cex=1.5)

  }
  dev.off()
  list(conv_ctcf=conv_ctcf, enr=enr)
}

fig_s_qc <- function(height_per_panel=85, width=520, karyo_zlim=c(-1, 1), karyo_cutoff=1, min_peak_loc=15.5)
{
  cells = rownames(sch_batch[order(sch_batch$batch), ])

  repli_score = rowSums(sch_tor_marg_n[cells, c('tor3', 'tor4')])
  repli_score = repli_score / min(repli_score, na.rm=T)
  vals = data.frame(coverage=rowSums(sch_glob_decay[cells,])/2,
    p_trans=sch_glob_decay[cells, 'ctrans']/rowSums(sch_glob_decay[cells,]),
    no_dig=sch_glob_decay[cells, 'c0']/rowSums(sch_glob_decay[cells,]),
    repli_score=repli_score,
    p_dup=rowMeans(sch_fend_dup[cells,]))

  hlines=list(coverage=c(sch_min_cov, sch_max_cov), p_trans=sch_max_trans, no_dig=sch_max_non_digest, repli_score=NULL, p_dup=NULL)


  fig_png(sprintf("%s/%s/s1_qcs.png",sch_fig_dir, rel_odir), width=width, height=(ncol(vals) + 1) * height_per_panel)
  layout(matrix(1:(ncol(vals)+1), ncol(vals)+1, 1))
  par(mar=c(0.5, 4, 0.5, 0.5))

  for (s in colnames(vals)) {
    cvals = as.vector(vals[,s])
    names(cvals) = cells

    .fig_s_cells_by_batch(cvals, s, sprintf("%s/%s/s1_qc_%s.png",sch_fig_dir, rel_odir, s), hlines=hlines[[s]], new_plot=F, plot_axis=F)
  }

  batch_counts = table(sch_batch[cells, 'batch'])
  batch_cumsum = cumsum(batch_counts)

  plot(seq_along(cells) / length(cells), col=NA, xlab="", xaxt='n', main="", ylab="labels")
  text(x=(c(0, batch_cumsum[-length(batch_cumsum)]) + batch_cumsum)/2, y=1, labels=names(batch_counts), adj=c(1, 0.5), srt=90, cex=1.2)
  dev.off()

  # karyo chrom z-score
  x = sch_calc_chrom_marg_z()
  karyo = x[["enr"]][cells, paste0("chr", sch_chroms)]
  bad_karyo_cells = rownames(karyo)[apply(abs(karyo), 1, max) > karyo_cutoff]

  fig_png(sprintf("%s/%s/s1_qc_karyo_all.png",sch_fig_dir, rel_odir), width=nrow(karyo)+length(bad_karyo_cells) + 50, height=600)
  layout(matrix(1:2, 1, 2), widths=c(nrow(karyo), length(bad_karyo_cells)))
  par(mar=c(0.5,4,0.5,0.5))
  image(pmin(pmax(karyo, karyo_zlim[1]), karyo_zlim[2]), col=colorRampPalette(diver_colspec)(256), zlim=karyo_zlim, main="", yaxt='n', xaxt='n')
  axis(2, at=seq(0, 1, length=ncol(sch_chrom_marg_z)), labels=gsub("chr", "", colnames(karyo)), las=2)
  batch_counts = table(sch_batch[rownames(karyo), 'batch'])
  batch_cumsum = cumsum(batch_counts)
  abline(v=c(0, batch_cumsum)/nrow(karyo), lwd=2)

  par(mar=c(0.5,1,0.5,0.5))
  image(pmin(pmax(karyo[bad_karyo_cells,], karyo_zlim[1]), karyo_zlim[2]), col=colorRampPalette(diver_colspec)(256), zlim=karyo_zlim, main="", yaxt='n', xaxt='n')
  dev.off()

  plot_legend(sprintf("%s/%s/s1_qc_karyo_leg.png",sch_fig_dir, rel_odir), karyo_zlim, colorRampPalette(diver_colspec)(256), "enr")

  dn = sch_glob_decay_res[sch_good_cells, ]

  dn_n = dn / rowSums(dn[,-1])
  dn_n[,1] = dn[,1]
  cdn = cor(dn_n)


  fig_png(sprintf("%s/%s/s1_qc_nodig_trans_decay_cor.png",sch_fig_dir, rel_odir), width=220, height=220)
  par(mar=c(4,4,1,1))
  diag(cdn) = NA
  plot(cdn[nrow(cdn),], type='l', col='darkgreen', lwd=2, ylim=c(-0.5, 0.5), main="", xlab="", ylab="pearson", xaxt='n')
  axis(1, c(1, seq(9, nrow(cdn)-1, by=8),nrow(cdn)), c("No dig", sapply(as.numeric(gsub("X", "", colnames(dn)[seq(9, ncol(dn)-1, by=8)])), n2str, 2), "Trans"), las=2, cex.axis=0.8)
  lines(cdn[1,], lwd=2, col='red')
  grid(col='black')
  legend("topleft", legend=c("Trans", "NoDig"), lty=1, lwd=2, col=c("darkgreen", "red"), cex=0.7, ncol=2, seg.len=1)
  dev.off()

  #near_f_vlines = c(sch_high_near_s_start, sch_high_near_mid_s)
  near_f_vlines = c(sch_near_s_start, sch_near_mid_s)

  fig_png(sprintf("%s/%s/s1_qc_decay_metric_density.png",sch_fig_dir, rel_odir), width=80*4, height=140)
  layout(matrix(1:5, 1, 5))
  par(mar=c(4,2,0.5,0.5))

  plot(density(sch_decay_metrics$mitotic_bin), main="", xlab="%M contacts", ylab="", lwd=1)
  grid(col='darkgray', lwd=0.5)
  abline(v=sch_mitotic_cutoff, lty=2, col='red')
  plot(density(sch_decay_metrics$high_near_f), main="", xlab="%near contacts", ylab="", lwd=1)
  grid(col='darkgray', lwd=0.5)
  abline(v=near_f_vlines, lty=2, col='red')
  plot(density(sch_decay_metrics$far_tightness), main="", xlab="far tightness", ylab="", xlim=c(0.02, 0.1), lwd=1)
  grid(col='darkgray', lwd=0.5)
  plot(density(sch_decay_metrics$early_f), main="", xlab="early cov", ylab="", lwd=1)
  grid(col='darkgray', lwd=0.5)
  plot(density(sch_decay_metrics$glob_max), main="", xlab="glob wm", ylab="", lwd=1)
  grid(col='darkgray', lwd=0.5)
  abline(v=min_peak_loc, lty=2, col='red')
  dev.off()

  scat_w = fig_factor * 420
  scat_h = fig_factor * 240
  alpha = 0.7

  fig1_decay_scatters(yfield="far_tightness", ylim=NULL, col_by_cond=F, vlines=near_f_vlines, width=scat_w, height=scat_h, alpha=alpha, cex=fig_factor * 0.3, vlwd=3 * fig_factor)
  fig1_decay_scatters(xfield="f_near_band", yfield="far_mu", ylim=NULL, col_by_cond=F, vlines=near_f_vlines, width=scat_w, height=scat_h, alpha=alpha, cex=fig_factor * 0.3, vlwd=3 * fig_factor)
  fig1_decay_scatters(xfield="f_near_band", y=sch_decay_metrics$early_f / min(sch_decay_metrics$early_f), yfield="early_f", ylim=NULL, col_by_cond=F, vlines=near_f_vlines, width=scat_w, height=scat_h, alpha=alpha, cex=fig_factor * 0.3, vlwd=3 * fig_factor)

  fig_png(sprintf("%s/%s/s1_qc_karyo_enr_per_chrom.png",sch_fig_dir, rel_odir), width=480, height=280)
  par(mfrow=c(4,5))
  par(mar=c(2,1,1,1))

  for (i in sch_chroms) {
    chrom = paste0('chr', i)
    hist(sch_chrom_marg_enr[, chrom], 500, xlim=c(-2,2), main=chrom)
    abline(v=c(-karyo_cutoff, karyo_cutoff), col='red', lty=2)
  }
  dev.off()

  # trans switcher
  if (!is.null(sch_trans_switchers)) {
    hits = 1:50
    
    y = log10(sch_trans_switchers$count+1)
    names(y) = sch_trans_switchers$n
    n_potential = sum(sch_trans_switchers$count)
    n_contact = sum(sch_trans_switchers$count * sch_trans_switchers$n)
    binom_exp = log10(n_potential * dbinom(x=hits, size=n_contact, prob=1/n_potential) + 1)

    fig_png(sprintf("%s/%s/s1_qc_trans_dups_barplot.png",sch_fig_dir, rel_odir), width=240, height=220)
    par(mar=c(4,4,1,1))
    plot(hits, y[hits], type='l', main="", xlab="#cells per contact", ylab="count (log10)")
    points(hits, y[hits], pch=19, col='red', cex=0.8)
    
    hits = hits[binom_exp > 0.1]
    lines(hits, binom_exp[hits], col='darkgray')
    points(hits, binom_exp[hits], pch=19, col='darkgray', cex=0.8)

    grid(col='black', lwd=0.5)
    
    legend("topright", legend=c('Obs', 'Exp'), pch=19, col=c('red', 'darkgray'), bty='n')
    
    dev.off()
  }

  ## nx = ceiling(sqrt(length(bs)))
  ## ny = ceiling(length(bs) / nx)

  ## fig_png(sprintf("%s/%s/s1_qc_batch_reads_per_contact.png",sch_fig_dir, rel_odir), width=120*nx, height=120*ny)
  ## layout(matrix(1:(nx*ny), nx, ny, byrow=T))
  ## par(mar=c(2,2,1,1))
  ## first=T
  ## for (b in bs) {
  ##   cb = colSums(sch_cis_contact_mul[ nmb[[b]], ])
  ##   tb = colSums(sch_trans_contact_mul[ nmb[[b]], ])
  ##   maxi = min(3, max(c(which(cb > 0), which(tb > 0))))

  ##   probs = data.frame(cis=cb/sum(cb), trans=tb/sum(tb))
  ##   barplot(t(probs[1:maxi, ]), beside=T, legend.text=first, col=c('red', 'darkgray'), cex.names=0.8, args.legend=list(bty='n', ncol=2), xlab="reads per contact", ylab="prob", main=b)
  ##   first=F
  ##   grid(col='black', lwd=0.25)
  ## }
  ## dev.off()

  # batch metrics table
  reads = read.table(sprintf("%s/reads_per_cell.txt", sch_table_dir))
  b = sch_batch
  b$good = is.element(rownames(b), sch_good_cells)
  b$reads = reads[rownames(b), ]

  write.table(as.data.frame(summarize(group_by(b, batch, batch_id, cond), n=length(good), n_good=sum(good), mrpc=mean(reads)/1e6)), sprintf("%s/batch_metrics.txt", sch_table_dir), quote=F, sep="\t")

  dec = sch_glob_decay_res[, grep("trans", colnames(sch_glob_decay_res), invert=T)]
  dec_dists = as.numeric(gsub("X", "", colnames(dec)))
  dind = dec_dists > sch_remove_near_dist

  dec = dec[, dind]
  dec_dists = dec_dists[dind]
  
  dec_n = dec / rowSums(dec)

  fig_png(sprintf("%s/%s/s1_which_max_decay_bin.png", sch_fig_dir,  rel_odir), width=220, height=280)
  h = hist(log2(dec_dists[apply(dec_n, 1, which.max)]), 200, col=darkblue_col, border=darkblue_col, main="", xlab="max bin")
  box()
  abline(v=log2(sch_glob_min_wm_dist_bin), lty=2, lwd=1*fig_factor)
  grid(col='black', lwd=0.5*fig_factor)
  dev.off()
  
}

fig_s_saturation <- function(min_reads=1e5, min_mols=5000, min_p_single=0, cex=0.5, col=darkblue_col, alpha=0.3, size=220, cex_lab=1.3)
{
  xa = as.matrix(sch_trans_contact_mul + sch_cis_contact_mul) / 2

  mols = xa %*% (1:ncol(xa))

  #p_single = xa[,1] / rowSums(xa)
  single_chains = read.table(sprintf("%s/p_singleton_chains.txt", sch_table_dir))
  sc = data.frame(p_single=single_chains[rownames(mols), ], batch=sch_batch[rownames(mols), 'batch'], col=unlist(sch_cond_colors[sch_batch[rownames(mols), 'cond']]))
  rownames(sc) = rownames(mols)
  
  reads = read.table(sprintf("%s/reads_per_cell.txt", sch_table_dir))
  #rownames(reads) = reads$V1
  reads = reads[rownames(mols), 'n']

  reads_lim = c(log10(min_reads), log10(max(reads)))
  mols_lim = c(log10(min_mols), log10(max(mols)))
  p_single_lim = c(min_p_single, max(sc$p_single))

  fig_png(sprintf("%s/%s/fig_s_reads_saturation1.png",sch_fig_dir, rel_odir), size, size)
  par(mar=c(4,4,1,1))
  plot(log10(reads), log10(mols), pch=19, cex=cex, col=addalpha(col, alpha), xlim=reads_lim, ylim=mols_lim, xlab="reads (log10)", ylab="mols (log10)", main="", cex.lab=cex_lab)
  abline(a=0, b=1, lwd=2, col='black', lty=3)
  grid(col='black', lwd=0.5)
  dev.off()

  #plot(log10(reads), p_single, pch=19, cex=cex, col=col, xlim=reads_lim, ylim=p_single_lim, xlab="reads (log10)", ylab="%single", main="", cex.lab=cex_lab)
  #grid(col='black', lwd=0.5)

  fig_png(sprintf("%s/%s/fig_s_reads_saturation2.png",sch_fig_dir, rel_odir), size, size)
  par(mar=c(4,4,1,1))
  plot(log10(reads), log2(reads/mols), pch=19, cex=cex, col=addalpha(col, alpha), xlim=reads_lim, xlab="reads (log10)", ylab="reads per mol (log2)", main="", cex.lab=cex_lab)
  grid(col='black', lwd=0.5)

  dev.off()

  # %singleton chains per batch
  fig_png(sprintf("%s/%s/fig_s_p_singleton_chains.png",sch_fig_dir, rel_odir), width=400, height=220)
  par(mar=c(4,4,1,1))

  bs = which(table(sch_batch$batch) >= 1)
  cells = rownames(sch_batch)[sch_batch$batch %in% names(bs)]
  sc = sc[cells, ]
  ccols = unique(sc[, c('batch', 'col')])

  boxplot(p_single ~ batch, sc, notch=F, col=as.character(ccols$col), xlab='batch', ylab='% singleton chains', pch=19, cex=0.5, ylim=c(0, 0.65))
  grid(col='black', lwd=0.5)

  dev.off()




}

.get_col_by_batch <- function(no_cond_col='darkgray')
{

  cols = rep(no_cond_col, nrow(sch_batch))
  names(cols) = rownames(sch_batch)

  cind = sch_batch$cond %in% names(sch_batch_colors)
  cols[cind] = unlist(sch_batch_colors[sch_batch[cind, 'cond']])

  cols
}

.fig_s_cells_by_batch <- function(vals, ylab, ofn, width=500, height=250, cex=0.2, ylim_q=c(0.002, 0.998), no_cond_col='black', hlines=NULL, new_plot=T, plot_axis=T)
{
  cols = .get_col_by_batch(no_cond_col=no_cond_col)
  cells = intersect(names(cols), names(vals))

  batch_counts = table(sch_batch[cells, 'batch'])
  batch_cumsum = cumsum(batch_counts)

  if (new_plot) {
    fig_png(ofn, width=width, height=height)
    par(mar=c(7, 4, 0.5, 0.5))
  }

  ylim = quantile(vals, ylim_q, na.rm=T)
  plot(pmin(pmax(vals[cells], ylim[1]), ylim[2]), pch=19, col=cols[cells], cex=cex, ylim=ylim, ylab=ylab, xlab="", xaxt='n', main="")
  abline(v=c(1, batch_cumsum), col='black')
  if (plot_axis) {
    axis(1, at=(c(0, batch_cumsum[-length(batch_cumsum)]) + batch_cumsum)/2, labels=names(batch_counts), las=2, cex.axis=0.8, tick=F)
  }

  if (!is.null(hlines)) {
    abline(h=hlines, lty=2)
  }

  if (new_plot) {
    dev.off()
  }

}


fig_s_subsample_cv_phasing <- function(mark_margin=0.1, n_partitions=10)
{
  cells = list()
  for (i in 1:5) {
    cells[[i]] = fig_s_2fold_cv_phasing(mark_margin=mark_margin, by_chrom=F, n_partitions=n_partitions, partitions=list(p1=seq(1, i), p2=seq(6, 5+i)))
  }
  cells
}

fig_s_2fold_cv_phasing <- function(mark_margin=0.1, by_chrom=T, partitions=list(p1=1:5, p2=6:10), n_partitions=10, chroms1=NULL, chroms2=NULL)
{
  chroms = gintervals(sch_chroms)
  chroms = chroms[order(chroms$end), ]

  if (by_chrom) {
    if (is.null(chroms1)) {
      chroms1 = as.character(chroms[seq(1, nrow(chroms), by=2), 'chrom'])
    }
    if (is.null(chroms2)) {
      chroms2 = as.character(chroms[seq(2, nrow(chroms), by=2), 'chrom'])
    }
    
    c1 = group_cells_by_metrics(plot_basic=F, plot_extra=F, criteria=sch_phasing_criteria, cells=rownames(sch_decay_metrics), chroms=chroms1)
    c2 = group_cells_by_metrics(plot_basic=F, plot_extra=F, criteria=sch_phasing_criteria, cells=rownames(sch_decay_metrics), chroms=chroms2)
  }
  else {
    c1 = group_cells_by_metrics(plot_basic=F, plot_extra=F, criteria=sch_phasing_criteria, cells=rownames(sch_decay_metrics), n_partitions=n_partitions, use_partitions=partitions[["p1"]])
    c2 = group_cells_by_metrics(plot_basic=F, plot_extra=F, criteria=sch_phasing_criteria, cells=rownames(sch_decay_metrics), n_partitions=n_partitions, use_partitions=partitions[["p2"]])
    
  }
  
  cells1 = c1$ord
  names(cells1) = rownames(c1)

  cells2 = c2$ord
  names(cells2) = rownames(c2)

  use_cells = intersect(names(cells1), names(cells2))
  

  if (by_chrom) {
    xlab = paste(gsub("chr", "", chroms1), sep="", collapse="-")
    ylab = paste(gsub("chr", "", chroms2), sep="", collapse="-")
  }
  else {
    xlab = paste0(partitions[["p1"]], collapse="-")
    ylab = paste0(partitions[["p2"]], collapse="-")
  }

  fig_png(sprintf("%s/%s/fig_s_2fold_phasing_CV_%s_%s_vs_%s.png",sch_fig_dir, rel_odir, ifelse(by_chrom, "chrom_halves", "contact_halves"), xlab, ylab), 320, 320)
  par(mar=c(4,4,2,1))

  plot(cells1[use_cells], cells2[use_cells], col="black", pch=19, cex=0.4, xlab=xlab, ylab=ylab, main=paste(length(use_cells), "cells"))
  first  = pmin(cells1[use_cells], cells2[use_cells])
  second = pmax(cells1[use_cells], cells2[use_cells])
  outliers = pmin(second - first, abs(second - first - length(use_cells))) >= mark_margin* length(use_cells)
  message(sprintf("max diff: %f\tout of %d cells\n", max(pmin(second - first, abs(second - first - length(use_cells)))/length(use_cells)), length(use_cells)))
  points(cells1[use_cells[outliers]], cells2[use_cells[outliers]], pch=19, cex=0.4, col='red')

  abline(v=cumsum(table(c1$group)), col='black', lty=1, lwd=1)
  abline(h=cumsum(table(c1$group)), col='black', lty=1, lwd=1)
  abline(a=mark_margin*length(use_cells), b=1, col='darkgray', lty=2, lwd=3)
  abline(a=-mark_margin*length(use_cells), b=1, col='darkgray', lty=2, lwd=3)

  dev.off()

  print(sum(outliers))

  use_cells
}

fig_s_chrom_metrics_phasing_by_cc <- function()
{

  res = chrom_phasing_stats_by_cc(odir=paste(sch_fig_dir, "paper", sep="/"), pointsize=pointsize, fres=fres)
}


fig_s_pool_knn_norm_cis_maps <- function(chroms=sch_chroms)
{
  commands = paste0("plot_knn_norm_cis_map(chr=\"", chroms, "\")", collapse=",")
  res = eval(parse(text=paste("gcluster.run(",commands,")")))

  shaman_colspec = c(rev(unlist(options("shaman.score_pal_neg_col"))), unlist(options("shaman.score_pal_0_col")), unlist(options("shaman.score_pal_pos_col")))

  shaman_breaks = c(-1 * rev(unlist(options("shaman.score_pal_neg_breaks"))), 0, unlist(options("shaman.score_pal_pos_breaks")))

  crp = colorRampPalette(shaman_colspec)(c(length(shaman_breaks)-1)*1000)
  plot_legend(ofn=sprintf("%s/%s/shaman_score_legend.png", sch_fig_dir, rel_odir), zlim=range(shaman_breaks), crp=crp[vals_to_cols(seq(min(shaman_breaks), max(shaman_breaks), length=1000), shaman_breaks, 1000)], "score")
}

plot_knn_norm_cis_map <- function(chr=19, coords=NULL, score_tn=paste0(pool_tn, "_lowerres_score"), bp_per_px=100e3, conv_ctcf=NULL, ctcf_chip_q=0.9, plot_loops=F, loop_size=3, loop_stroke=1, plot_also_domains=F, extend_loop_by=1e4, min_d=60, points=NULL, band=NULL, point_size=0.1, grid_every=NULL)
{
  library(shaman)
  message(packageVersion("shaman"))
  message(packageVersion("misha"))
  .update_shaman_options()
  
  if (is.null(coords)) {
    cint = gintervals(chr)
  }
  else {
    cint = gintervals(chr, coords[1], coords[2])
  }

  if (is.null(band)) {
    band = c(cint$start - cint$end, 0)
  }

  if (plot_loops & is.null(conv_ctcf)) {
    message("getting ctcf sites")
    conv_ctcf = .get_conv_ctcf_with_tor(pool_q_cutoff=0, ctcf_chip_q=ctcf_chip_q, pool_loop_extend=extend_loop_by, min_d_score=min_d)
  }
  curr_loops = conv_ctcf[conv_ctcf$chrom1 == paste0("chr", chr) & conv_ctcf$start2 > conv_ctcf$start1, ]

  width = round((cint$end - cint$start) / bp_per_px)

  if (is.null(points)) {
    message("getting points")
    points = gextract(score_tn, intervals=gintervals.2d(chr), colnames="score", band=band)
  }

  fig_fn=sprintf("%s/%s/fig_s_chr%s%s%s.png",sch_fig_dir, rel_odir, as.character(chr), ifelse(is.null(coords), "", sprintf("_%s-%s", n2str(coords[1], 2), n2str(coords[2], 2))), ifelse(plot_loops, "_ctcf_loops", ""))

  message("plotting")

  if (plot_loops | !is.null(grid_every) | plot_also_domains) {
    map_score = shaman_gplot_map_score(points_score=points, interval_range=cint, add_axis=T, point_size=point_size)

    if (plot_loops) {
      loops = data.frame(x=(curr_loops$start1 + curr_loops$start2)/2, y=(curr_loops$start2 - curr_loops$start1)/2)
      map_score = map_score + geom_point(data=loops, aes(x, y), size=loop_size, stroke=loop_stroke, col='black', shape=1)
    }
    
    if (plot_also_domains) {
      d = .sch_get_pool_domains()
      d = d[d$chrom == paste0('chr', chr) & d$end >= cint$start & d$start <= cint$end, ]
      ds = data.frame(x=c(d$start, (d$end + d$start)/2), y=c(rep(0, nrow(d)), (d$end - d$start)/2), xend=c((d$end + d$start)/2, d$end), yend=c((d$end - d$start)/2, rep(0, nrow(d))))
      map_score = map_score + geom_segment(data=ds, mapping=aes(x=x, y=y, xend=xend, yend=yend), color='black', size=1, linetype=1)
    }

    if (!is.null(grid_every)) {
      xs = seq(cint$start, cint$end, by=grid_every)
      gridl = data.frame(x=c(xs, xs), y=rep(0, 2*length(xs)), xend=c((xs + cint$end)/2, c(xs + cint$start)/2), yend=c((cint$end - xs)/2, (xs - cint$start)/2))
      map_score = map_score + geom_segment(data=gridl, mapping=aes(x=x, y=y, xend=xend, yend=yend), color='black', size=0.2, linetype=1)
    }
    fig_png(fig_fn, width=width, height=width)

    #map_score + ggsave(fig_fn, width=width, height=(5/6) * width, units='px')
    print(map_score)
    dev.off()
  }
  else {
    shaman_plot_map_score_with_annotations(genome="mm9", points_score=points, interval_range=cint, point_size=point_size, add_genes=F, add_ideogram=F, fig_fn=fig_fn, fig_width=width, fig_height=(5/6) * width)
  }

  list(points=points, conv_ctcf=conv_ctcf)
}


fig_s_group_decay_and_ensemble_ref <- function(width=320, height=240, ref_nm="hic.ES.tg_esc", chroms=sch_chroms, rdec_n=NULL, hg19_mit_dec_n=NULL, hg19_mit_nm="hic.K562.mitotic", m_d=2**c(21, 23.5), n_mitotic=8, lwd=2*fig_factor)
{

  m = sch_decay_metrics[sch_decay_metrics$valid, ]
  dec = sch_glob_decay_res[rownames(m), -ncol(sch_glob_decay_res)]

  dec_dists = as.numeric(gsub("X", "", colnames(dec)))
  dec = dec[, dec_dists > sch_remove_near_dist]
  dec_dists = dec_dists[ dec_dists > sch_remove_near_dist]

  dec_n = dec / rowSums(dec)

  if (is.null(rdec_n)) {
    scope = gintervals(chroms)
    rdec = gcis_decay(ref_nm, sch_cis_breaks, scope , scope)
    rdec = rdec[,1][attr(rdec, 'breaks')[-1] > sch_remove_near_dist - 1]
    rdec_n = rdec / sum(rdec)
  }

  if (is.null(hg19_mit_dec_n)) {
    gdb.init("/net/mraid14/export/tgdata/db/tgdb/hg19/trackdb", rescan=T)
    scope = gintervals(c(1:22, 'X'))

    hg19_dec = gcis_decay(hg19_mit_nm, sch_cis_breaks, scope , scope)
    hg19_dec = hg19_dec[,1][attr(hg19_dec, 'breaks')[-1] > sch_remove_near_dist - 1]
    hg19_mit_dec_n = hg19_dec / sum(hg19_dec)
    gdb.init("/net/mraid14/export/tgdata/db/tgdb/mm9/trackdb", rescan=T)
  }

  fig_png(sprintf("%s/%s/fig_s_decay_ens_ref_mit_naumova.png",sch_fig_dir, rel_odir), width=width, height=height)

  f_mit    = rowSums(dec_n[, dec_dists > m_d[1] & dec_dists <= m_d[2]])
  mit_cells = tail(names(f_mit)[order(f_mit)], n=n_mitotic)

  par(mar=c(5,3,0.5,0.5))

  plot(colMeans(dec_n), type='l', lwd=lwd, ylim=c(0,0.033), xlab="", ylab="", xaxt='n', col='black')
  li = seq(1, ncol(dec_n), by=5)
  axis(1, at=li, labels=sapply(dec_dists, n2str, 1)[li], cex=0.6, las=2)

  lines(rdec_n, lwd=lwd, col='blue', lty=3)
  lines(hg19_mit_dec_n, lwd=lwd, col='red', lty=3)
  lines(colMeans(dec_n[mit_cells,]), lwd=lwd, col='orange')

  grid(col='black', lwd=0.5)
  legend("topleft", legend=c("Pool", "Ens mES", "K562 M", "Top M"), lty=c(1,3,3,1), lwd=lwd, col=c('black', 'blue', 'red', 'orange'), bty='n', cex=1*fig_factor)

  dev.off()


  fig_png(sprintf("%s/%s/fig_s_decay_by_phases.png",sch_fig_dir, rel_odir), width=width, height=height)
  par(mar=c(3,3,0.5,0.5))

  plot(colMeans(dec_n), type='l', lwd=lwd, ylim=c(0,0.03), xlab="", ylab="", xaxt='n', col='darkgray')
  li = seq(1, ncol(dec_n), by=6)
  axis(1, at=li, labels=sapply(dec_dists, n2str, 1)[li], cex=0.7, las=1)

  lines(rdec_n, lwd=lwd, col='black')
  lines(hg19_mit_dec_n, lwd=2, col='orange', lty=3)
  lines(colMeans(dec_n[m$group == 1,]), lwd=lwd/2, col='orange')
  lines(colMeans(dec_n[m$group == 2,]), lwd=lwd/2, col='blue')
  lines(colMeans(dec_n[m$group == 3,]), lwd=lwd/2, col='red')
  lines(colMeans(dec_n[m$group == 4,]), lwd=lwd/2, col='darkgreen')
  grid(col='black', lwd=0.5)

  legend("topleft", legend=c("Pool", "Ens mES", "K562 M"), lty=c(1,1,3), lwd=lwd, col=c('darkgray', 'black', 'orange'), bty='n', cex=0.9*fig_factor)

  legend("topright", legend=paste(c("M (", "G1 (", "Early-S (", "Late-S/G2 ("), table(m$group), ")", sep=""), lwd=lwd/2, col=c('orange', 'blue', 'red', 'darkgreen'), bty='n', cex=0.9*fig_factor)

  dev.off()

  list(rdec_n=rdec_n, hg19_mit_dec_n=hg19_mit_dec_n)
}

fig_s_domains_trans_cluster <- function(km=NULL, vfunc="weighted.sum", rna_tns=c("rna.129_Cast.ES.pf_rep1", "rna.129_Cast.ES.pf_rep2"), acol='red', bcol='darkgray')
{

  for (i in seq_along(rna_tns)) {
    gvtrack.create(paste0("rna_gpm", i), rna_tns[i], "global.percentile.max")
  }

  if (is.null(km)) {
    km = cluster_pool_hic_domains_by_inter_domain_contacts(fn_pref=sprintf("%s/%s/fig_s_domains_trans_cluster_%s",sch_fig_dir, rel_odir, vfunc), vfunc=vfunc)
  }

  km$d$clus = ifelse(km$d$cluster == 2, "A", "B")
  print(table(km$d$clus))

  km$d$rna = gextract(sprintf("-log2(1-pmax(%s))", paste0("rna_gpm", seq_along(rna_tns), collapse=",")), intervals=km$d, iterator=km$d, colnames="x")$x


  fig_png(sprintf("%s/%s/fig_s_domains_AB_stats_boxplot_%s.png",sch_fig_dir, rel_odir, vfunc), 540, 180)
  par(mfrow=c(1, 5))
  par(mar=c(2, 2, 2, 0.5))
  #boxplot(len ~ clus, km$d, main="length", notch=T, col=c(acol, bcol))
  boxplot(tor ~ clus, km$d, main="ToR", notch=T, col=c(acol, bcol))
  boxplot(k4me3 ~ clus, km$d, main="H3K4me3", notch=T, col=c(acol, bcol), ylim=c(0, 0.12))
  boxplot(k27me3 ~ clus, km$d, main="H3K27me3", notch=T, col=c(acol, bcol), ylim=c(0, 2.2))
  boxplot(laminB1 ~ clus, km$d, main="LaminB1", notch=T, col=c(acol, bcol))
  
  boxplot(rna ~ clus, km$d, main="RNA-Seq", notch=T, col=c(acol, bcol), ylim=c(2, 21))
  dev.off()

  fig_png(sprintf("%s/%s/fig_s_AB_dom_lengths.png",sch_fig_dir, rel_odir), 220, 180)
  al = density(km$d$len[km$d$clus == 'A'])
  bl = density(km$d$len[km$d$clus == 'B'])
  par(mar=c(2,2,1,1))

  plot(al$x, al$y, type='l', lwd=2, col=acol, main="", xlim=range(c(al$x, bl$x)), ylim=range(c(al$y, bl$y)), xlab='length', ylab='density')
  lines(bl$x, bl$y, lwd=2, col=bcol)
  grid(col='black', lwd=0.5)
  legend("topright", legend=c('A', 'B'), lty=1, lwd=2, col=c(acol, bcol), bty='n')
  dev.off()
  km
}

fig_s_decay_kmeans <- function(k=12, plot_mean_cmap=T, mean_cmap_rel_dir="paper/cmaps/decay_clust_log2_by_0.125", mean_cmap_scope=gintervals.2d(6, 0, -1, 6, 0, -1), name="_log2_by_0.125", zlim=c(0, 0.025), remove_no_digest=T, remove_trans=T, width=2400, height=500, colspec=seq_colspec, ncols=256, glob_dens_lim=c(1e-18, 5e-16), cmap_colspec=c(seq_colspec, 'black'))
{
  decay = sch_glob_decay_res[sch_good_cells, ]
  if (remove_no_digest) {
    decay = decay[, -1]
  }
  if (remove_trans) {
    decay = decay[, -ncol(decay)]
  }

  col_dists = as.numeric(gsub("X", "", colnames(decay)))
  decay = decay[, col_dists > sch_remove_near_dist]
  col_dists = col_dists[col_dists > sch_remove_near_dist]

  decay = decay / rowSums(decay)

  km_f = TGLKMeans_wrapper(decay, "/tmp/clus.tab", k)

  #hc = hclust(dist(km_f$centers))
  #km = order(hc$order)[km_f$cluster]
  #km = order(rowSums(km_f$centers[, col_dists >= 2e6 & col_dists <= 50e6]), decreasing=T)[km_f$cluster]
  #corder = order(c(12, 8, 4, 1, 10, 7, 9, 3, 6, 11, 2, 5))
  #km = corder[km_f$cluster]

  km_f$repli = rowSums(sch_tor_marg_n[names(km_f$cluster), c(3,4)])

  fig_png(sprintf("%s/%s/fig_s_decay_clusts%s_repli_boxplot.png",sch_fig_dir, rel_odir, name), width=540, height=320)
  boxplot(repli ~ cluster, data.frame(repli=km_f$repli, cluster=km_f$cluster), ylab="repli-score", notch=T, pch=19, col='darkgray', cex=0.5)
  grid(col='black', lwd=0.5)
  dev.off()

  km = km_f$cluster
  print(table(km))

  #width = length(km)
  fig_png(sprintf("%s/%s/fig_s_decay_clusts%s.png",sch_fig_dir, rel_odir, name), width=width, height=height)
  layout(matrix(1:2, 2, 1), heights=c(1, 9))

  
  par(mar=c(0,5,1,1))
  m = sch_decay_metrics
  m$clus = km[rownames(m)]
  image(as.matrix(as.numeric(m[ order(m$clus + 1e-5 * as.numeric(m$group)), 'group'])), col=mini_qualitative_colors[1:5], main="", xaxt='n', yaxt='n')
  clus_lines = cumsum(table(km)) / length(km)
  abline(v=clus_lines, lwd=2, col='black')
  
  par(mar=c(5,5,1,1))

  image(pmin(pmax(as.matrix(decay[order(km),]), zlim[1]), zlim[2]), col=colorRampPalette(colspec)(ncols), zlim=zlim, xaxt='n', yaxt='n')
  ax_at = seq(0, 1, len=ncol(decay))
  ax_c = seq(1, ncol(decay), by=4)
  ax_labs = sapply(col_dists, n2str, 2)
  if (!remove_trans) {
    ax_labs = c(ax_labs, "trans")
  }

  abline(v=clus_lines, lwd=2, col='black')

  axis(2, ax_at[ax_c], ax_labs[ax_c], las=2, cex.axis=1.4)
  axis(1, (c(0, clus_lines[-length(clus_lines)]) + clus_lines)/2, 1:k, cex.axis=1.4, lwd=0)

  dev.off()

  system(sprintf("mkdir -p %s/%s",sch_fig_dir, mean_cmap_rel_dir))
  clus = data.frame(x=km)
  rownames(clus) = names(km_f$cluster)
  plot_clus_decays_and_stats(only_decays=T, mean_cmap_rel_dir, clus, other_col='darkgrey', clus_col='red', width=fig_factor*130, height=fig_factor*100)
  plot_clus_decays_and_stats(only_scatters=T, mean_cmap_rel_dir, clus, other_col='darkgrey', clus_col='red', width=fig_factor*180, height=fig_factor*100)

  if (plot_mean_cmap) {
    sch_plot_mean_mat(sprintf("%s/%s",sch_fig_dir, rel_odir, mean_cmap_rel_dir), clus, dens_white_quantile=0.02, scope=mean_cmap_scope, fn_pref="decay_c", glob_dens_lim=glob_dens_lim, colspec=cmap_colspec)

  }

  km
}


fig_s_insu <- function(width=380, height=150, sw=101, cex=0.4, collapse_q=c(0.005, 0.995), col_by_cond=T, insu_data=NULL, dixon_insu_data=NULL, rebuild=F, filter_cov_domains=F, haploid=F)
{
  mv = sch_decay_metrics[ sch_decay_metrics$valid, ]

  if (is.null(dixon_insu_data)) {
    insu_name = sprintf("_Dixon_tad_borders_gt%s%s", n2str(sch_remove_near_dist), ifelse(filter_cov_domains, "_good_cov", ""))
    dixon_insu_data = insu_by_cc_grouped(sch_decay_metrics, insu_name=insu_name)
  }

  if (is.null(insu_data)) {
    insu_name = sprintf("_tad_borders_gt%s%s", n2str(sch_remove_near_dist), ifelse(filter_cov_domains, "_good_cov", ""))
    insu_data = insu_by_cc_grouped(sch_decay_metrics, insu_name=insu_name)
  }

  dy = -dixon_insu_data$mean_ins_per_cell
  y = -insu_data$mean_ins_per_cell

  cols = get_cond_cols(haploid)

  fig_png(sprintf("%s/%s/fig_s_dixon_bord_insu_sw%d%s.png",sch_fig_dir, rel_odir, sw, ifelse(filter_cov_domains, "_covD", "")), width=width, height=height)
  par(mar=c(3,3,0.5,0.5))
  plot(dy, cex=cex, pch=19, col=cols, main="", xlab="", ylab="insulation")
  grid(col='black', lwd=0.5)
  lines(cyclic_rollmean(dy, sw), lwd=2)
  lines(cyclic_rollmean(y, sw), lwd=2, col='red')

  legend("topright", legend=c("Dixon", "Ours"), lty=1, lwd=2, col=c('black','red'), bty='n', cex=0.7)
  dev.off()

  fig_png(sprintf("%s/%s/fig_s_dixon_vs_our_bord_scatter%s.png",sch_fig_dir, rel_odir, ifelse(filter_cov_domains, "_covD", "")), width=height, height=height)
  par(mar=c(3,3,0.5, 0.5))
  cols = darkblue_col
  plot(y, dy, pch=19, cex=0.2, col=addalpha(cols, 0.4), main="", xlab="ours", ylab="Dixon", asp=1)
  grid(col='black', lwd=0.5)
  abline(a=0, b=1, lwd=1)
  dev.off()


  fig_png(sprintf("%s/%s/fig_s_mean_insu_vs_near.png", sch_fig_dir, rel_odir), width=height, height=height)
  par(mar=c(3,3,0.5, 0.5))
  cols = darkblue_col
  plot(mv[names(y), 'f_near_band'], dy, pch=19, cex=0.4*fig_factor, col=addalpha(cols, 0.2), main="", xlab="%near", ylab="mean insu")
  grid(col='black', lwd=0.5)
  dev.off()

  list(insu_data=insu_data, dixon_insu_data=dixon_insu_data)
}


fig_s_hap_vs_dip_rna_seq_over_genes <- function()
{
    options(gmax.data.size=1e9)

    spec_params_fn_hap = "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/serum_2ic/serum_2i_271116_params.r"
    sch_load_config(spec_params_fn_hap)
    fn = sprintf("%s/%s/fig_s_dip_hap_rna_seq.png", sch_fig_dir, rel_odir)

    genes = gintervals.load("intervs.global.tss")
    genes = genes[genes$chrom %in% paste0("chr", sch_chroms), c("chrom","start","end")]

    gvtrack.create("hap_rna_1", sch_rna_tn_rep1, "global.percentile.max")
    gvtrack.create("hap_rna_2", sch_rna_tn_rep2, "global.percentile.max")
    hap_rna = gextract("-log2(1-((hap_rna_1 + hap_rna_2)/2))", intervals=genes, iterator=genes, colnames="rna")

    spec_params_fn_dip = "/net/mraid14/export/data/users/lubling/datasets/scell/results/esh/hyb_mm9_2i_is/hyb_mm9_2i_is_300317_params.r"
    sch_load_config(spec_params_fn_dip)
    gvtrack.create("dip_rna_1", sch_rna_tn_rep1, "global.percentile.max")
    gvtrack.create("dip_rna_2", sch_rna_tn_rep2, "global.percentile.max")
    dip_rna = gextract("-log2(1-((dip_rna_1 + dip_rna_2)/2))", intervals=genes, iterator=genes, colnames="rna")

    lmfit = lm(hap_rna$rna ~ dip_rna$rna)
    cex_val = 2.5
    png(fn, width=900, height=900)
    par(mar=c(10,10,10,10))
    plot(panel.first=grid(col="darkgray", lty="solid"),
      dip_rna$rna, hap_rna$rna, col=addalpha(brewer.pal(11, "RdYlBu")[11],0.2),
      log="xy",
      main=sprintf("Haploid vs. diploid mouse ES cell RNA-seq (R^2=%s)",round(summary(lmfit)$r.squared,3)),
      xaxt="n", yaxt="n", xlab="", ylab="", pch=19,
      cex.axis=cex_val, cex.lab=cex_val, cex.main=cex_val)
    # abline(lmfit)
    axis(side=1, cex.axis=cex_val, mgp=c(3,2,0))
    axis(side=2, cex.axis=cex_val, mgp=c(3,2,0))
    title(xlab="Diploid RNA-seq", cex.lab=cex_val, line=6)
    title(ylab="Haploid RNA-seq", cex.lab=cex_val, line=5)

    dev.off()

}

fig_s_trans_ab <- function(sw=101, width=380, height=180, cex=0.4, haploid=F, collapse_q=c(0.005, 0.995))
{
  mv = sch_decay_metrics[ sch_decay_metrics$valid, ]

  cis_ab = calc_ab_enrich(sch_ab_counts[rownames(mv),  grep("cis", colnames(sch_ab_counts))], "cis")
  trans_ab = calc_ab_enrich(sch_ab_counts[rownames(mv),  grep("trans", colnames(sch_ab_counts))], "trans")

  ylim = quantile(c(trans_ab$deplAB, cis_ab$deplAB), collapse_q, na.rm=T)
 
 
  fig_png(sprintf("%s/%s/fig_s_trans_ab_comp_sw%d.png",sch_fig_dir, rel_odir, sw), width=width, height=height)
  par(mar=c(3,3,0.5,0.5))

  cols = get_cond_cols(haploid)
  cols [ trans_ab$deplAB < ylim[1] | trans_ab$deplAB > ylim[2]] = 'red'
  plot(trans_ab$deplAB, cex=cex, pch=19, col=cols, main="", xlab="", ylab="AB comp", ylim=ylim)
  grid(col='black', lwd=0.5)
  lines(cyclic_rollmean(trans_ab$deplAB, sw), lwd=2)
  lines(cyclic_rollmean(cis_ab$deplAB, sw), lwd=2, lty=3)

  legend("topleft", legend=c("Trans", "Cis"), lty=c(1,3), lwd=2, col='black', bty='n', cex=0.7)
  dev.off()


}

fig_s_mitotic_ins <- function(mi="scell.nextera.pool_good_hyb_2i_idx_sort_es_group_1_post_m_ins_discard8192_300000s", msi="scell.nextera.pool_good_hyb_2i_idx_sort_es_group_1_post_m_shuffle_ins_discard8192_300000s", chrom="chr1", xlim = c(32e6, 38e6), ylim=c(1, 3.6), chroms=sch_chroms, scale_distrib=F, height=160)
{
  ins = gextract(ins_tn, mi, msi, intervals=gintervals(chroms), iterator=2e4, colnames=c('gi', 'mi', 'msi'))

  fig_png(sprintf("%s/%s/fig_s_mitotic_ins_example.png",sch_fig_dir, rel_odir), 280, height)
  par(mar=c(2,2,0.5,0.5))
  cins = ins[as.character(ins$chrom) == chrom, ]
  plot (cins$start, -cins$gi, type='l', col='black', xlim=xlim, ylim=ylim, lwd=2)
  lines(cins$start, -cins$mi, col='darkorange1', lwd=2)
  lines(cins$start, -cins$msi, col='red', lwd=2)
  grid(col='black', lwd=0.5)
  legend("topleft", legend=c("M", "M (shuf)", "Pool"), col=c('darkorange1', 'red', 'black'), lwd=2, lty=1, bty='n', ncol=3)
  dev.off()

  fig_png(sprintf("%s/%s/fig_s_mitotic_ins_distrib_%s.png",sch_fig_dir, rel_odir, ifelse(scale_distrib, "scaled", "raw")), 100, height)
  par(mar=c(2,2,0.5,0.5))
  if (scale_distrib) {
    plot(density(scale(-ins$msi, scale=F), na.rm=T), col='red', xlim=scale(ylim, scale=F), lwd=2, main="", xlab="insulation")
    lines(density(scale(-ins$mi, scale=F), na.rm=T), col='darkorange1', lwd=2)
    lines(density(scale(-ins$gi, scale=F), na.rm=T), col='black', lwd=2)
  }
  else {
    plot(density(-ins$msi, na.rm=T), col='red', xlim=ylim, lwd=2, main="", xlab="insulation")
    lines(density(-ins$mi, na.rm=T), col='darkorange1', lwd=2)
    lines(density(-ins$gi, na.rm=T), col='black', lwd=2)
  }

  grid(col='black', lwd=0.5)
  dev.off()

  ins

}

fig_s_trans_examples <- function(odir=sprintf("%s/%s/cmaps/trans_examples",sch_fig_dir, rel_odir), bp_per_px=1e6, cex=0.2, alpha=1, grid_every_bp=0, pool_vfunc="weighted.sum", pool_min_contacts=2, pool_binsize=5e5, chroms=matrix(c(8, 16, 10, 11, 2, 12), 3, 2, byrow=T))
{
  system(paste("mkdir -p", odir))
  
  nms = paste0("scell.nextera.1CD", c('X4_242', 'U_354', 'X2_402', 'X2_95', 'X4_224', 'X4_62', 'X1_284', 'U_428', 'X2_104', 'S2_148', 'X4_234', 'X4_353', 'X4_453', 'U_255', 'U_484', 'S1_178', 'X3_203', 'X4_232', 'X4_72'))

  for (i in 1:nrow(chroms)) {
    print(i)
    sch_plot_genome_wide(pool_tn, sprintf("%s/pool_%d_vs_%d_%s_min%d_bs%s.png", odir, chroms[i, 1], chroms[i, 2], pool_vfunc, pool_min_contacts, n2str(pool_binsize, 0)), chroms=chroms[i,], binsize=pool_binsize, vfunc=pool_vfunc, min_contacts=pool_min_contacts, legend_ofn=sprintf("%s/pool_%d_vs_%d_%s_min%d_bs%s_leg.png", odir, chroms[i, 1], chroms[i, 2], pool_vfunc, pool_min_contacts, n2str(pool_binsize, 0)))

    for (nm in nms) {
      message(nm)
      plot_trans_helper(nm, chroms[i, ], odir=odir, bp_per_px=bp_per_px, cex=cex, alpha=alpha, grid_every_bp=grid_every_bp)
    }
  }
}

fig_s_decay_and_repli_scatters <- function(domain_d=2**c(17, 19), near_d=sch_near_dists, m_d=sch_mitotic_dists, far_d=sch_far_dists, panel_size=180, w2h=2, cex=0.5, alpha=0.5, col='blue', pch=19, cex_lab=1.3, log_data=T, mitotic_to_far_d=2**c(17, 25.625),
                                           spec_nms=c(
                                             'scell.nextera.1CDX4_242',
                                             'scell.nextera.1CDX1_186',
                                             'scell.nextera.1CDS2_51',
                                             'scell.nextera.1CDX1_391',
                                             'scell.nextera.1CDU_555',
                                             'scell.nextera.1CDS2_543'), mark_top_n_mit=8, show_spec_nms=F, mitotic_like_cutoff=0.19)
{

  mv = sch_decay_metrics[sch_good_cells, ]

  dn = sch_glob_decay_res[rownames(mv), ]
  dn = dn[, -ncol(dn)]

  dn_dists = as.numeric(gsub("X", "", colnames(dn)))
  vind = dn_dists >= sch_remove_near_dist

  dn = dn[, vind]
  dn_dists = dn_dists[vind]

  dn_n = dn / rowSums(dn)

  df = data.frame(
    domain = rowSums(dn_n[, dn_dists > domain_d[1] & dn_dists <= domain_d[2]]),
    mit    = rowSums(dn_n[, dn_dists > m_d[1] & dn_dists <= m_d[2]]),
    near   = rowSums(dn_n[, dn_dists > near_d[1] & dn_dists <= near_d[2]]),
    far    = rowSums(dn_n[, dn_dists > far_d[1] & dn_dists <= far_d[2]]),
    m_to_far    = rowSums(dn_n[, dn_dists > mitotic_to_far_d[1] & dn_dists <= mitotic_to_far_d[2]]))

  rownames(df) = rownames(mv)

  x = df$mit
  y = df$near

  if (log_data) {
    df = log2(df)
  }

  df$repli  = mv$early_f / min(mv$early_f)

  cols = rep('#999999', nrow(df))
  names(cols) = rownames(mv)
  strong_m_cutoff = quantile(df$m_to_far, 5/6)
  #cols[df$m_to_far >= strong_m_cutoff] = 'blue'

  cols[df$mit >= ifelse(log_data, log2(mitotic_like_cutoff), mitotic_like_cutoff)] = '#333333'
  cols[tail(rownames(df)[order(df$mit)], n=mark_top_n_mit)] = 'red'

  print(table(cols))

  mit_lab = sprintf("mit %s-%s", n2str(m_d[1], 1), n2str(m_d[2], 1))
  dom_lab = sprintf("domain %s-%s", n2str(domain_d[1], 1), n2str(domain_d[2], 1))
  near_lab = sprintf("near %s-%s", n2str(near_d[1], 1), n2str(near_d[2], 1))
  far_lab = sprintf("far %s-%s", n2str(far_d[1], 1), n2str(far_d[2], 1))
  mit_to_far = sprintf("M-far %s-%s", n2str(mitotic_to_far_d[1], 1), n2str(mitotic_to_far_d[2], 1))
  if (!is.null(spec_nms)) {
    nm_c = sapply(spec_nms, function(nm) { which(rownames(df) == nm) })
  }


  .scatter_figs_s_helper = function(v1, s1, l1, v2, s2, l2) {
    fig_png(sprintf("%s/%s/fig_s_decay_%s_vs_%s_scatter%s%s.png",sch_fig_dir, rel_odir, s1, s2, ifelse(show_spec_nms, "_spec", ""),  ifelse(log_data, "_log", "")), panel_size*w2h, panel_size)
    par(mar=c(4,4,1,1))
    plot(v1, v2, pch=pch, cex=cex, col=addalpha(cols, alpha), xlab=l1, ylab=l2, cex.lab=cex_lab, main="")
    ind = cols == 'blue'
    points(v1[ind], v2[ind], pch=pch, cex=cex, col=addalpha('blue', alpha))

    grid(col='black', lwd=0.5)
    if (show_spec_nms) {
      text(v1[nm_c], v2[nm_c], seq_along(nm_c), pos=3, col='magenta', cex=2)
    }
    dev.off()
  }
  .scatter_figs_s_helper(df$mit, "mit", mit_lab, df$domain, "domain", dom_lab)
  .scatter_figs_s_helper(df$mit, "mit", mit_lab, df$near, "near", near_lab)
  .scatter_figs_s_helper(df$mit, "mit", mit_lab, df$far, "far", far_lab)
  .scatter_figs_s_helper(df$mit,    "mit",    mit_lab, df$repli, "repli", "repli-score")
  .scatter_figs_s_helper(df$domain, "domain", dom_lab, df$repli, "repli", "repli-score")
  .scatter_figs_s_helper(df$far,    "far",    far_lab, df$repli, "repli", "repli-score")
  #.scatter_figs_s_helper(df$repli, "repli", "repli-score", df$mit,    "mit",    mit_lab)
  #.scatter_figs_s_helper(df$repli, "repli", "repli-score", df$domain, "domain", dom_lab)
  #.scatter_figs_s_helper(df$repli, "repli", "repli-score", df$far,    "far",    far_lab)

  fig_png(sprintf("%s/%s/fig_s_repli_ee_vs_ll_marg_cov.png",sch_fig_dir, rel_odir), 330, 330)
  par(mar=c(3,3,1,1))

  mv_p = mv[ permute(1:nrow(mv)), ]
  plot(sch_tor_marg_n[rownames(mv_p), 4], sch_tor_marg_n[rownames(mv_p), 1], pch=19, cex=cex, xlab="EE", ylab="LL", asp=1, col=addalpha(mini_qualitative_colors[as.numeric(mv_p$group)], alpha))
  ind = as.numeric(mv$group) == 2
#  points(sch_tor_marg_n[ind, 4], sch_tor_marg_n[ind, 1], pch=19, cex=cex, col=addalpha(mini_qualitative_colors[2], alpha))
                                        #plot(sch_tor_marg_n[rownames(mv), 4], sch_tor_marg_n[rownames(mv), 1], pch=19, cex=cex, xlab="EE", ylab="LL", asp=1, col=cols)
  #cind = cols != '#999999'
  #points(sch_tor_marg_n[rownames(mv)[cind], 4], sch_tor_marg_n[rownames(mv)[cind], 1], pch=19, cex=0.4, col=cols[cind])
  grid(col='black', lwd=0.5)
  dev.off()

  # Fig 1E distribution of the mitotic band strength on all cells, indicating a threshold for mitotic-like behavior
  width = 190
  height= 270
  fig_png(sprintf("%s/%s/fig_1e_mitotic_band_distrib.png",sch_fig_dir, rel_odir), width, height)
  dd = density(2**df$mit)
  plot(dd, main="", xlab="%mitotic", ylab="density", lwd=2)
  lines(dd$x[dd$x > mitotic_like_cutoff], dd$y[dd$x > mitotic_like_cutoff], col='palevioletred1', lwd=2)
  #abline(v=mitotic_like_cutoff, col='blue', lty=2, lwd=2)
  grid(col='darkgray')
  dev.off()

  # Fig 1F distribution of repli-score
  fig_png(sprintf("%s/%s/fig_1f_repli_score_distrib.png",sch_fig_dir, rel_odir), width, height)
  repli_score = 2**df$repli
  repli_score = repli_score / min(repli_score)
  plot(density(repli_score), main="", xlab="repli-score", ylab="density", lwd=3)
  grid(col='darkgray')
  dev.off()

  df$cols = cols
  df
}

#####
fig_s_chrom_chrom_total_contacts <- function(reg=1e-3, psize=100, norm_cutoff=0.3, xchrs = c("chr1", "chr6", "chr18"), ychrs = c("chr2", "chr7", "chr19")

)
{
  x = sch_trans_pairs_contact_count[sch_good_cells, grep("chrM|chrY", colnames(sch_trans_pairs_contact_count), invert=T)
    ]
  x_n = x / rowSums(x)

  x_nn = t(log2(t(x_n + reg) / colMeans(x_n + reg)))

  chr1 = gsub(" .*", "", colnames(x), perl=T)
  chr2 = gsub(".* ", "", colnames(x), perl=T)

  chroms = gintervals.all()
  clen = chroms$end / 10^6
  names(clen) = chroms$chrom

  x_d = 1e3 * t(t(x) / (clen[chr1] * clen[chr2]))

  fig_png(sprintf("%s/%s/fig_s_trans_pair_examples.png",sch_fig_dir, rel_odir), psize*length(xchrs), psize*length(ychrs))


  layout(matrix(1:(length(xchrs) * length(ychrs)), length(xchrs), length(ychrs)))
  par(mar=c(2,2,1,1))

  for (i in seq_along(xchrs)) {
    for (j in seq_along(ychrs)) {
      plot(density(log10(x_d[, paste(xchrs[i], ychrs[j])])), main="", xlab="", ylab="", lwd=2)
      abline(v=norm_cutoff, col='red', lty=3, lwd=2)
      #legend("topleft", legend=gsub("chr", "", paste(ychrs[j], xchrs[i], sep=':')), col=NA, bty='n')
      #if (j == length(ychrs)) {
      #  title(xlab=xchrs[i], cex=0.7)
      #}
      #if (i == 1) {
      #  title(ylab=ychrs[j], cex=0.7)
      #}
      grid(col='black', lwd=0.5)
    }
  }
  dev.off()

  fig_png(sprintf("%s/%s/fig_s_trans_pairs_interact_hm.png",sch_fig_dir, rel_odir), 280, 280)
  par(mar=c(2,2,1,1))
  m = matrix(0, length(sch_chroms), length(sch_chroms))
  rownames(m) = colnames(m) = paste0("chr", sch_chroms)
  m[cbind(chr1, chr2)] = colMeans(x_d >= 10^norm_cutoff)

  image(m, col=colorRampPalette(seq_colspec)(256), zlim=c(0, 1), xaxt='n', yaxt='n')
  axis(1, at=seq(0, 1, len=nrow(m)), sch_chroms, lwd=0, cex.axis=0.7)
  axis(2, at=seq(0, 1, len=nrow(m)), sch_chroms, lwd=0, cex.axis=0.7)

  dev.off()

  plot_legend(sprintf("%s/%s/fig_s_trans_pairs_interact_hm_leg.png",sch_fig_dir, rel_odir), c(0, 1), colorRampPalette(seq_colspec)(256), "frac")

}

fig_s_repli <- function(colspec=diver_colspec, zlim=c(-0.3, 0.3), rebuild=F)
{
  for (s in c('tor', 'trans_A_score')) {
    dd = domain_marg_cov_plots(order_by=s, ord_only=F, group_size=0, odir=sprintf("%s/%s",sch_fig_dir, rel_odir), colspec=colspec, name=paste0(colspec, collapse="-"), zlim=zlim, rebuild=rebuild)
  }

  #dd = domain_marg_cov_plots(ord_only_by_tor=F, group_size=0, odir=sprintf("%s/%s",sch_fig_dir, rel_odir), colspec=colspec, name=paste0(colspec, collapse="-"), zlim=zlim, order_by="none", clus_method="hclust")

}

compile_features_table <- function(col_names=c('ord', 'group', 'f_mitotic_band', 'f_near_band', 'far_mu'), insu_data=NULL, loops_d=NULL, group_names=c('post-M', 'G1', 'early-S', 'late-S/G2', 'pre-M'), calc_loops=T, calc_insu=T, haploid=F)
{
  m = sch_batch
  nms = rownames(m)

  m$cell_nm = gsub("scell.nextera.", "", nms)
  
  m$f_no_digestion = sch_glob_decay[nms, 'c0']/rowSums(sch_glob_decay[nms,])
  m$mean_f_fend_dup = rowMeans(sch_fend_dup[nms, ])
  m$max_chrom_cov_abber = apply(abs(sch_chrom_marg_enr[nms, ]), 1, max)
  
  align = read.table(sprintf("%s/trans_contacts_corr.tab", sch_table_dir), header=T, stringsAsFactors=F)
  rownames(align) = align$nm
  align = align[nms,]

  m$f_trans = align$ntrans / align$n
  m$trans_align = align$c_align
  
  m$total_contacts = round(rowSums(sch_chrom_stat[nms,])/2)

  cis_ab = calc_ab_enrich(sch_ab_counts[nms,  grep("cis", colnames(sch_ab_counts))], "cis")
  trans_ab = calc_ab_enrich(sch_ab_counts[nms,  grep("trans", colnames(sch_ab_counts))], "trans")
  m$cis_ab_comp   = cis_ab$deplAB
  m$trans_ab_comp = trans_ab$deplAB

  m$passed_qc = ifelse(is.element(nms, sch_good_cells), 1, 0)

  # only for cells passing QC
  if (calc_insu) {
    if (is.null(insu_data)) {
      insu_data = insu_by_cc_grouped(sch_decay_metrics, insu_name=sprintf("_tad_borders_gt%s", n2str(sch_remove_near_dist)))
    }
    
    m[names(insu_data$mean_ins_per_cell), 'mean_insu'] = -insu_data$mean_ins_per_cell
  }

  if (calc_loops) {
    if (is.null(loops_d)) {
      if (haploid) {
        loops_d = fig3_conv_ctcf_loops_by_cc_haploid()
      }
      else {
        loops_d = fig3_conv_ctcf_loops_by_cc()
      }
    }
    l_en = loops_d$x_enr

    m[rownames(l_en), 'loop_enr_EE'] = l_en$EE_EE
    m[rownames(l_en), 'loop_enr_LL'] = l_en$LL_LL
  }


  for (s in col_names) {
    m[rownames(sch_decay_metrics), s] = sch_decay_metrics[, s]
  }
  m[rownames(sch_decay_metrics), 'repli_score'] = sch_decay_metrics$early_f / min(sch_decay_metrics$early_f)

  m$group = group_names[as.numeric(m$group)]

  write.table(m, sprintf("%s/%s_features_table.txt", sch_table_dir, ifelse(haploid, "haploids", "diploids")), quote=F, sep="\t", row.names=F)
  
  
  list(insu_data=insu_data, loops_d=loops_d)
  
}


fig_s_fends_pairs_orientation <- function(stats=NULL, mindist=100, maxdist=1e4, log_step=0.125, cols=list(f00='#e41a1c', f01='#377eb8', f10='#4daf4a', f11='#984ea3'), leg=list(f00='right', f01='conv', f10='diver', f11='left'), lwd=2*fig_factor)
{
  if (is.null(stats)) {
    stats = fends_pair_orientation(mindist=mindist, maxdist=maxdist, log_step=log_step)
  }

  fig_png(sprintf("%s/%s/fig_s_fends_distrib_by_orientation_%s_to_%s_by_%f.png", sch_fig_dir, rel_odir, n2str(mindist, 0), n2str(maxdist, 0), log_step), 280, 280)

  stats = stats[, -ncol(stats)] # remove last partial bin
  
  d = as.numeric(gsub("b", "", colnames(stats)))
  
  plot(d, stats[1, ], col=NA, ylim=range(stats), main="", xlab="distance (log2)", ylab="#contacts")
  
  apply(cbind(unlist(cols[rownames(stats)]), stats), 1, function(x) { lines(d, x[-1], col=x[1], lwd=lwd) })
  legend("bottomright", legend=unlist(leg[rownames(stats)]), col=unlist(cols[rownames(stats)]), lty=1, lwd=lwd, bty='n')
  grid(col='black', lwd=0.5*fig_factor)
  dev.off()

  stats
}

collect_fendchain_stats <- function(cis_cutoff=1e3, fns_regexp="1CDU|1CDX", redb_dir="/net/mraid14/export/tgdata/db/tgdb/mm9/seq/redb", re="GATC", rebuild=F, only_good=F)
{
  ifns = unlist(lapply(sch_base_dir, list.files, pattern="^fendchain$", recursive=T, full.names=T))
  ifns = ifns[ grep(fns_regexp, ifns, perl=T)]
  cell_nms = paste0("scell.nextera.", gsub(".", "_", sapply(gsub("/fendchain", "", ifns), basename), fixed=T))

  if (only_good) {
    ind = is.element(cell_nms, sch_good_cells)
    ifns = ifns[ind]
    cell_nms = cell_nms[ind]
  }
  ofns = paste0(sch_rdata_dir, "/fc_report_", cell_nms)
  
  if (rebuild) {
    commands = paste0("system(\"perl ", Sys.getenv("PIPELINE_HOME"), "/map3c/TG3C/fendchain_report.pl ", " ", redb_dir, " ", re, " ", ifns, " ", cis_cutoff, " ", ofns, "\")", collapse=",")

    res = eval(parse(text=paste("gcluster.run(",commands,")")))
  }

  counts = do.call("rbind", lapply(ofns, function(ifn) { x = read.delim(ifn, header=T, row.names=NULL); x = as.data.frame(arrange(summarize(group_by(x, n_trans, n_far), count=sum(count)), n_trans, n_far)); x$count = x$count / sum(x$count); x }))
#      bins_cov = as.data.frame(arrange(summarize(group_by(x, chrom1, start1, end1), total=sum(count)), chrom1, start1))
  counts_n =  as.data.frame(arrange(summarize(group_by(counts, n_trans, n_far), freq=mean(count)), n_trans, n_far))
  m = matrix(0, 5, 5, dimnames=list(paste0("cis", 0:4), paste0("trans", 0:4)))
  m[ cbind(paste0("cis", counts_n$n_far), paste0("trans", counts_n$n_trans))]  = counts_n$freq
  write.table(m, sprintf("%s/%s/fendchain_%s_cells_cis_gt%s_trans_report.txt", sch_fig_dir, rel_odir, ifelse(only_good, "good", "all"), n2str(cis_cutoff, 0)), quote=F, sep="\t")
}

# compare insulation of different tracks
fig_s_compare_tracks_insu <- function(tns=c(
                                        "scell.nextera.pool_good_hyb_2i_idx_sort_es",
                                        "scell.nextera.pool_good_hyb_2i_es",
                                        "scell.nextera.pool_good_hyb_idx_sort_es",
                                        "scell.nextera.pool_good_hyb_serum_es", 
                                        "scell.nextera.pool_good_es_c",
                                        "scell.nextera.pool_good_serum_es",
                                        "scell.nextera.pool_good_2i_es",
                                        "hic.ES.tg_esc"),
                                      nms=c('dip', 'dip_sort', 'dip_nosort', 'dip_serum', 'hap', 'hap_serum', 'hap_2i', 'ens'), pairs=matrix(c(1, 8, 1, 5, 1, 4, 2, 3, 7, 6), 5, 2, byrow=T), chroms=sch_chroms, ins_res=1e3, sample_every=10e3, cex=0.1*fig_factor, alpha=0.4, lim=c(1,9), scale=ins_scale, discard_below=1000, size=220)
{
  ints = gintervals(chroms)
  ins_tns = c()
  for (tn in tns) {
    ins_tn = sprintf("%s_ins_%ds", tn, scale)
    ins_tns = c(ins_tns, ins_tn)
    if (!gtrack.exists(ins_tn)) {
      message("producing ins track ", ins_tn)
      gtrack.2d.gen_insu_track(tn, scale, res=ins_res, new_track=ins_tn, min_diag_d=discard_below)
    }
  }

  x = eval(parse(text=paste0("gextract(\"",paste0(ins_tns, collapse="\",\""),"\", intervals=ints, iterator=sample_every, colnames=nms)")))

  if (is.null(lim)) {
    lim = range(-x[, nms], na.rm=T)
  }

  for (i in 1:nrow(pairs)) {
    i1 = pairs[i,1]
    i2 = pairs[i,2]
    fig_png(sprintf("%s/%s/fig_s_insu_%s_vs_%s_s%s_r%s.png", sch_fig_dir, rel_odir, nms[i1], nms[i2], n2str(scale, 0), n2str(sample_every, 0)), size, size)
    par(mar=c(4,4,1,1))
    plot(-x[, nms[i1]], -x[, nms[i2]], pch=19, cex=cex, main="", xlab=nms[i1], ylab=nms[i2], col=addalpha(darkblue_col, alpha), asp=1, xlim=lim, ylim=lim)
    points(-x[, nms[i1]], -x[, nms[i2]], pch=19, cex=cex, col=addalpha(darkblue_col, alpha/30))
    legend("topleft", legend=sprintf("cor=%.2f", cor(x[, nms[i1]], x[, nms[i2]], use='pair')), col=NA)
    #abline(h=mean(-x[, nms[i2]], na.rm=T), col='black', lty=2)
    #abline(v=mean(-x[, nms[i1]], na.rm=T), col='black', lty=2)

    #abline(h=-quantile(x[, nms[i2]], ins_dom_thresh, na.rm=T), col='black', lty=2)
    #abline(v=-quantile(x[, nms[i1]], ins_dom_thresh, na.rm=T), col='black', lty=2)
    #abline(a=0, b=1, col='darkgray')

    grid(col='black', lwd=0.5*fig_factor)
    dev.off()
  }

}

########################
# A-score plots
fig_s_a_score_ideograms <- function(tn=pool_tn, ref_tn="scell.nextera.pool_good_hyb_2i_idx_sort_es",  min_cis=2e6, vfunc="weighted.sum", min_total_d_contacts=20, colspec=rev(brewer.pal(n=11, "RdYlBu")[c(1:4, 8:11)]), d=NULL, nshades=1000, chroms=paste0("chr", sch_chroms), calc_frac_a=F, pch_col=darkblue_col, alpha=0.2, cex=0.5 * fig_factor, rebuild=F, lim=c(0.54, 0.78), black_domains=NULL) 
{
  stopifnot(!calc_frac_a)
  if (is.null(d)) {

    if (tn == ref_tn) {
      cis_d = calc_domains_a_score(tn=ref_tn, vfunc=vfunc, rebuild=rebuild, use_trans=F)
      trans_d = calc_domains_a_score(tn=ref_tn, vfunc=vfunc, rebuild=rebuild, use_trans=T)
    }
    else {
      ref_d = sch_get_ab_tagged_domains(raw_tn=tn, vfunc=vfunc, by_ref_trans_a_score_tn=ref_tn, rebuild=rebuild)
      cis_d = calc_domains_a_score(tn=tn, vfunc=vfunc, d=ref_d, rebuild=rebuild, use_trans=F)
      trans_d = calc_domains_a_score(tn=tn, vfunc=vfunc, d=ref_d, rebuild=rebuild, use_trans=T)
      
    }

    d = merge(cis_d, trans_d)
  }

  if (is.null(lim)) {
    if (calc_frac_a) {
      lim = c(0, 1) }
    else {
      lim = range(d[, c('trans_A_score', 'cis_A_score')], na.rm=T)
    }
  }

  type_lab = ifelse(calc_frac_a, "P(A)", "A-score")
  
  v = vals_to_cols(pmin(pmax(d$trans_A_score, lim[1]), lim[2]), lim, nshades)
  names(v) = paste(d$chrom, d$start, sep="_")
  cols = colorRampPalette(colspec)(nshades)
  if (!is.null(black_domains)) {
    v[black_domains] = nshades + 1
    cols = c(cols, 'black')
  }
  sch_plot_domain_score_idiogram(d, v, sprintf("%s/%s/trans_A_score_idiogram_%s.png", sch_fig_dir, rel_odir, ifelse(calc_frac_a, "pA", "A_score")), cols, "", chrs=chroms, title=paste(tn, type_lab, "trans"))
  
  v = vals_to_cols(d$cis_A_score, lim, nshades)
  names(v) = paste(d$chrom, d$start, sep="_")
  cols = colorRampPalette(colspec)(nshades)
  if (!is.null(black_domains)) {
    v[black_domains] = nshades + 1
    cols = c(cols, 'black')
  }

  sch_plot_domain_score_idiogram(d, v, sprintf("%s/%s/cis_A_score_idiogram_%s.png", sch_fig_dir, rel_odir, ifelse(calc_frac_a, "pA", "A_score")), cols, "", chrs=chroms, title=paste(tn, type_lab, "cis"))

  plot_legend(sprintf("%s/%s/cis_tran_pA_idiograms_%s_leg.png", sch_fig_dir, rel_odir, ifelse(calc_frac_a, 'pA', 'A_score')), lim, colorRampPalette(colspec)(nshades), type_lab)


  #orig A/B by trans contacts idiogram
  orig_ab_col = grep("orig_ab", colnames(d))
  stopifnot(length(orig_ab_col) == 1)
  sch_plot_domain_score_idiogram(d, ifelse(d[, orig_ab_col] == 'A', 1, 2), sprintf("%s/%s/trans_AB_clustered_idiogram.png", sch_fig_dir, rel_odir), c(diver_colspec[11], diver_colspec[1]), chrs=chroms, legend_nms=c('A', 'B'), title=sprintf("%s (%d domains)", pool_tn, nrow(d)))

  valid = d$A + d$B >= min_total_d_contacts
  fig_png(sprintf("%s/%s/cis_trans_A_score_scatter_%s.png", sch_fig_dir, rel_odir, ifelse(calc_frac_a, 'pA', 'A_score')), 340, 340)
  plot(d$trans_A_score[valid], d$cis_A_score[valid], pch=19, cex=cex, xlab=paste("trans", type_lab), ylab=paste("cis", type_lab), col=addalpha(pch_col, alpha), asp=1)
  grid(col='black', lwd=0.5*fig_factor)
  dev.off()

  fig_png(sprintf("%s/%s/cis_trans_A_score_densities_%s.png", sch_fig_dir, rel_odir, ifelse(calc_frac_a, 'pA', 'A_score')), 340, 340)
  plot(density(d$trans_A_score[valid], na.rm=T), col='red', lwd=2, xlab=type_lab, main="")
  lines(density(d$cis_A_score[valid], na.rm=T), col='blue', lwd=2)
  grid(col='black', lwd=0.5*fig_factor)
  legend("topleft", legend=c('trans', 'cis'), col=c('red', 'blue'), lwd=2, bty='n', cex=2*fig_factor)
  dev.off()
  
  d
}

###########
# compare domains A-score based on pool of singles domains
fig_s_compare_tracks_a_score <- function(ref_tn="scell.nextera.pool_good_hyb_2i_idx_sort_es", ref_nm='dip',
                                         tns=c(
                                           "scell.nextera.pool_good_hyb_2i_idx_sort_es",
                                           "scell.nextera.pool_good_hyb_2i_es",
                                           "scell.nextera.pool_good_hyb_idx_sort_es",
                                           "scell.nextera.pool_good_hyb_serum_es",
                                           "scell.nextera.pool_good_es_c",
                                           "scell.nextera.pool_good_serum_es",
                                           "scell.nextera.pool_good_2i_es",
                                           "hic.ES.tg_esc",
                                           "scell.nextera.pool_good_hyb_2i_idx_sort_es_group_2_g1",
                                           "scell.nextera.pool_good_hyb_2i_idx_sort_es_group_3_early_s",
                                           "scell.nextera.pool_good_hyb_2i_idx_sort_es_group_4_mid_s_g2"),
                                         nms=c('dip', 'dip_sort', 'dip_nosort', 'dip_serum', 'hap', 'hap_serum', 'hap_2i', 'ens', 'dip_2i_g1', 'dip_2i_early_s', 'dip_2i_mid_s_g2'), pairs=matrix(c(1, 8, 1, 5, 1,4, 2, 3, 7, 6, 9, 10, 9, 11, 10, 11), 8, 2, byrow=T), vfunc="weighted.sum", cex=0.2*fig_factor, alpha=0.5, lim=NULL, col_dict=list('A'='red', 'B'='black'), size=240, rebuild=F)
{
  x = calc_domains_a_score(tn=ref_tn, vfunc=vfunc, d=NULL, use_trans=T, rebuild=rebuild)
  
  for (i in seq_along(tns)) {
    message("calc A-score for ", nms[i])
    as = calc_domains_a_score(tn=tns[i], vfunc=vfunc, d=x[, c('chrom', 'start', 'end', 'ab')], use_trans=T, rebuild=rebuild)
    x[, nms[i]] = NA
    x[rownames(as), nms[i]] = as$trans_A_score
  }

  if (is.null(lim)) {
    lim = range(x[, nms], na.rm=T)
  }

  cols = unlist(col_dict[x$ab])
  
  for (i in 1:nrow(pairs)) {
    i1 = pairs[i,1]
    i2 = pairs[i,2]
    fig_png(sprintf("%s/%s/fig_s_a_score_by_%s_%s_vs_%s_%s.png", sch_fig_dir, rel_odir, ref_nm, nms[i1], nms[i2], vfunc), size, size)
    par(mar=c(4,4,3,3))
    plot(x[, nms[i1]], x[, nms[i2]], pch=19, cex=cex, main=sprintf("cor=%.2f", cor(x[, nms[i1]], x[, nms[i2]], use='pair')), xlab=nms[i1], ylab=nms[i2], col=addalpha(cols, alpha), asp=1, xlim=lim, ylim=lim)
    legend("topleft", legend=paste("Dip. ", names(col_dict)), col=unlist(col_dict), pch=19, bty='n')

    #abline(h=mean(-x[, nms[i2]], na.rm=T), col='black', lty=2)
    #abline(v=mean(-x[, nms[i1]], na.rm=T), col='black', lty=2)

    #abline(h=-quantile(x[, nms[i2]], ins_dom_thresh, na.rm=T), col='black', lty=2)
    #abline(v=-quantile(x[, nms[i1]], ins_dom_thresh, na.rm=T), col='black', lty=2)
    #abline(a=0, b=1, col='darkgray')

    grid(col='black', lwd=0.5*fig_factor)
    dev.off()
  }

}

# Haploid - Diploid pc contig table
fig_s_hap_dip_pc_contig_table <- function(size=280 * fig_factor, cex=0.5 * fig_factor, zlim=c(0, 0.02), colspec=c('white', seq_colspec))
{
  hap = read.table("~/schic_res/serum_2ic/tables/domain_ord_by_mean_g4_A_assoc_with_new_cells.txt", header=T)
  dip = read.table("~/schic_res/hyb_mm9_2i_is/tables/domain_ord_by_mean_g4_A_assoc.txt", header=T)
  nn = gintervals.neighbors.wrap(hap, dip, maxneighbors=100, maxdist=0)
  hdt = table(nn$pc_1, nn$pc_2) / nrow(nn)

  fig_png(sprintf("%s/%s/fig_s_hap_vs_dip_pc.png", sch_fig_dir, rel_odir), size, size)
  par(mar=c(4,4,1,1))
  image(pmin(hdt, zlim[2]), zlim=zlim, col=colorRampPalette(colspec)(256), xaxt='n', yaxt='n', main="", xlab="hap PC", ylab="dip PC")

  n_pc = max(dip$pc)

  axis(1, at=seq(0, 1, len=n_pc), labels=1:n_pc, tick=F)
  axis(2, at=seq(0, 1, len=n_pc), labels=1:n_pc, tick=F)
  
  #abline(v=0.5/n_pc+seq(0, 1, len=n_pc+1), lwd=0.5*fig_factor)
  #abline(h=0.5/n_pc+seq(0, 1, len=n_pc+1), lwd=0.5*fig_factor)

  dev.off()

  plot_legend(sprintf("%s/%s/fig_s_hap_vs_dip_pc_leg.png", sch_fig_dir, rel_odir), zlim, colorRampPalette(colspec)(256), 'fraction')
  

}

##
# params fit haploid. see below call for diploid
plot_domains_weighted_a_score_across_groups <- function(ref_tn="scell.nextera.pool_good_hyb_2i_idx_sort_es", min_cis_dist=2e6, ra=NULL, vfunc="weighted.sum", min_pool_dom_comp_contacts=20, min_domain_contacts=4, min_unique_domains_contacted=3, size=170, pch_col=darkblue_col, alpha=0.2, cex=0.4*fig_factor, min_obs_per_group=10, frac_tab_zmax=0.4, use_downsampled=T, freq_bin=0.03, ab_colspec=diver_dark_mid_colspec, freq_colspec=seq_colspec, freq_lim=c(0.54, 0.78), nshades=1000, order_doms_by="mean_g4", n_pc=20, n_clusters=30, cells_per_slice=101, ab_lim=NULL, ab_score_bin=0.04, mean_sd_by_slice=NULL, use_frac_c=F, sd_bp_ylim=c(0.01, 0.07), rebuild=F)
{
  if (is.null(ra)) {
    ra = calc_domains_weighted_a_score_across_cells(ref_tn=ref_tn, min_cis_dist=min_cis_dist, vfunc=vfunc, min_domain_contacts=min_domain_contacts, min_unique_domains_contacted=min_unique_domains_contacted, use_frac_c=use_frac_c, rebuild=rebuild)
  }

  d = ra$ref_d
  ind = !is.na(d$trans_A_score) & d$A + d$B >= min_pool_dom_comp_contacts

  d = d[ind, ]
  mean_mu = ra$mu_A[ind, ]
  total = ra$total_A[ind, ]
  n_d = ra$n_dom_A[ind, ]
  ds_mu = ra$ds_mu_A[ind, ]

  if (use_downsampled) {
    mu = ds_mu
  }
  else {
    mu = mean_mu
  }
  mu = mu[, rownames(sch_decay_metrics)]

  v = d$trans_A_score
  v = v - min(v)
  v = v / max(v)
  d$color = colorRampPalette(ab_colspec)(151)[1 + floor(v * 100)]

  
  ofn_pref = sprintf("%s/%s/fig_domains_a_score_across_groups_by_%s_%s_minPool%d_minTotal%d_minContDoms%d_%s", sch_fig_dir, rel_odir, pool_tn, vfunc, min_pool_dom_comp_contacts, min_domain_contacts, min_unique_domains_contacted, ifelse(use_downsampled, "_ds", "_mean"))


  fig_png(sprintf("%s_tor_vs_As.png", ofn_pref), 2*size, 2*size)
  plot(d[!is.na(d$tor), 'trans_A_score'], d[!is.na(d$tor), 'tor'], pch=19, cex=cex, col=addalpha(pch_col, alpha), xlab="A-score", ylab="mean ToR")
  grid(col='black', lwd=0.5*fig_factor)
  dev.off()
  
  fig_png(sprintf("%s_var_vs_dlen.png", ofn_pref), 2*size, 2*size)
  plot(d$end - d$start, apply(mu, 1, sd, na.rm=T), pch=19, cex=cex, col=addalpha(d$color, alpha), xlab="domain length", ylab="mean A-score sd")
  grid(col='black', lwd=0.5*fig_factor)
  dev.off()


  
  fig_png(sprintf("%s_mean_of_mean_vs_ds.png", ofn_pref), 2*size, 2*size)
  plot(apply(mean_mu, 1, mean, na.rm=T), apply(ds_mu, 1, mean, na.rm=T), pch=19, cex=cex, col=addalpha(d$color, alpha), xlab="mean meanA", ylab="mean dsA")
  grid(col='black', lwd=0.5*fig_factor)
  dev.off()

  fig_png(sprintf("%s_a_assoc_vs_a_score.png", ofn_pref), 2*size, 2*size)
  plot(d$trans_A_score, apply(mu, 1, mean, na.rm=T), pch=19, cex=cex, col=addalpha(pch_col, alpha), xlab="trans A-score", ylab="mean A-assoc")
  grid(col='black', lwd=0.5*fig_factor)
  dev.off()


  #ab_breaks = seq(ifelse(is.null(ab_lim), floor(min(mu, na.rm=T) / ab_score_bin) * ab_score_bin, ab_lim[1]), ifelse(is.null(ab_lim), ceiling(max(mu, na.rm=T) / ab_score_bin) * ab_score_bin, ab_lim[2]), by=ab_score_bin)


  valid = rowSums(!is.na(mu)) >= min_obs_per_group * 2

  fig_png(sprintf("%s_mu_sd_vs_mean.png", ofn_pref), 2*size, 2*size)
  mu_m = apply(mu, 1, mean, na.rm=T)
  mu_sd = apply(mu, 1, sd, na.rm=T)

  ab_breaks = seq(ifelse(is.null(ab_lim), quantile(mu_m[valid], 0.01, na.rm=T), ab_lim[1]), ifelse(is.null(ab_lim), quantile(mu_m[valid], 0.99, na.rm=T), ab_lim[2]), len=7)

  ab_cols = colorRampPalette(ab_colspec)(length(ab_breaks) - 1)

  ab_score_bin = ab_breaks[2] - ab_breaks[1]

  plot(mu_m[valid], mu_sd[valid], pch=19, cex=cex, col=ab_cols[cut(mu_m[valid], breaks=ab_breaks, labels=1:(length(ab_breaks)-1))], xlab="mean mu", ylab="sd mu", main=paste(sum(valid), "domains"))
  #points(mu_m[!valid], mu_sd[!valid], pch=19, cex=cex, col='black')
  grid(col='black', lwd=0.5*fig_factor)
  dev.off()

  mu_ind = split(1:sum(valid), cut(pmin(pmax(d[valid, 'trans_A_score'], min(ab_breaks)), max(ab_breaks)), breaks=ab_breaks, include.lowest=T))
  
  if (is.null(mean_sd_by_slice)) {
    mean_sd_by_slice = sapply(mu_ind, function(v) { colMeans(cyclic_rollapply(mu[valid, ][v, ], width=cells_per_slice, by=1, func_str="sd", na.rm=T), na.rm=T) })
  }

  fig_png(sprintf("%s_mu_sd_by_As_and_cc_b%.2f.png", ofn_pref, ab_score_bin), 2*size, 2*size)
  plot(mean_sd_by_slice[, 1], col=NA, ylim=range(mean_sd_by_slice, na.rm=T), main="", xlab="phased position", ylab="SD (mean)")
  abline(v=cumsum(table(sch_decay_metrics$group)))
  apply(rbind(ab_cols, mean_sd_by_slice), 2, function(v) { lines(v[-1], col=v[1], lwd=2) })
  grid(col='black', lwd=0.5 * fig_factor)
  legend("bottomleft", legend=ab_breaks[1:ncol(mean_sd_by_slice)], col=colorRampPalette(ab_colspec)(ncol(mean_sd_by_slice)), lwd=2, ncol=2, cex=0.8)  
  dev.off()

  gr_ind = split(1:ncol(mu), sch_decay_metrics[colnames(mu), 'group'])
  gr_ind = gr_ind[c('2', '3', '4')]
  
  muv = mu[valid, ]
  for (di in seq_along(mu_ind)) {
    sd_df = NULL
    for (gi in seq_along(gr_ind)) {
      sd_df = rbind(sd_df, data.frame(mean_sd=apply(muv[ mu_ind[[di]], gr_ind[[gi]] ], 1, sd, na.rm=T), group=names(gr_ind)[gi]))
    }
    mw23 = wilcox.test(sd_df[sd_df$group == '2', 'mean_sd'], sd_df[sd_df$group == '3', 'mean_sd'], alternative="greater")
    mw34 = wilcox.test(sd_df[sd_df$group == '3', 'mean_sd'], sd_df[sd_df$group == '4', 'mean_sd'], alternative="greater")
    
    message(sprintf("%s\t2-3: %g\t3-4: %g", names(mu_ind)[di], mw23$p.value, mw34$p.value))
      
    fig_png(sprintf("%s_sd_by_As_boxplot_s%d.png", ofn_pref, di), size, size*2)
    boxplot(mean_sd ~ group, sd_df, varwidth=T, notch=T, col=ab_cols[di], main=di, ylab="mean SD", outline=F, ylim=sd_bp_ylim)
    print(table(sd_df$group))
    dev.off()
  } 

  freq_breaks = seq(ifelse(is.null(freq_lim), floor(min(mu, na.rm=T) / freq_bin) * freq_bin, freq_lim[1]), ifelse(is.null(freq_lim), ceiling(max(mu, na.rm=T) / freq_bin) * freq_bin, freq_lim[2]), by=freq_bin)

  grps = split(sch_good_cells, sch_decay_metrics[sch_good_cells, 'group'])
  
  fig_png(sprintf("%s_phase_group_metrics.png", ofn_pref), size * length(grps),  size * 3)
  layout(matrix(1: (length(grps) * 3), 3, length(grps)))

  sd_lim = range(apply(mu, 1, sd, na.rm=T), na.rm=T)
  mean_lim = range(apply(mu, 1, mean, na.rm=T), na.rm=T)

  mu_freq = list()
  g_mu = list()
  g_sd = list()
  valid_doms = c()
  for (g in seq_along(grps)) {
    g_nm = names(grps)[g]
    plot(ecdf(rowSums(!is.na(mu[, grps[[g]]]))), do.points=F, xlab='#valid', ylab='ecdf', main=sprintf("group %s (%d)", g_nm, length(grps[[g]])), lwd=2*fig_factor)
    grid(col='black', lwd=0.5 * fig_factor)
    
    plot(ecdf(apply(total[, grps[[g]]], 1, median, na.rm=T)), do.points=F, xlab='#per d (med)', ylab='ecdf', main="", lwd=2*fig_factor)
    grid(col='black', lwd=0.5 * fig_factor)

    g_valid = rowSums(!is.na(mu[, grps[[g]]])) >= min(min_obs_per_group, length(grps[[g]]))
    plot(apply(mu[g_valid, grps[[g]]], 1, mean, na.rm=T), apply(mu[g_valid, grps[[g]]], 1, sd, na.rm=T), pch=19, cex=cex, xlab='mean sc-Ascore', ylab='sd sc-Ascore', main=sprintf("%d Doms", sum(g_valid)), xlim=mean_lim, ylim=sd_lim, col=addalpha(pch_col, alpha))
    grid(col='black', lwd=0.5 * fig_factor)

    g_mu[[g_nm]] = apply(mu[g_valid, grps[[g]]], 1, mean, na.rm=T)
    g_sd[[g_nm]] = apply(mu[g_valid, grps[[g]]], 1, sd, na.rm=T)
    
    mu_b = apply(mu[g_valid, grps[[g]]], 1, function(v) { table(cut(v, breaks=freq_breaks)) } )
    mu_freq[[g_nm]] = t(mu_b) / colSums(mu_b, na.rm=T)

    if (g_nm == "2") {
      valid_doms = rownames(mu)[g_valid]
    }
    else if (g_nm == "3" | g_nm == "4") {
      valid_doms = intersect(valid_doms, rownames(mu)[g_valid])
    }
      
  }
  
  dev.off()

  fig_png(sprintf("%s_comp_grps.png", ofn_pref), size * 3, size * 2)
  par(mar=c(4,4,1,1))
  layout(matrix(1:6, 2, 3))
  gind = matrix(as.character(c(2,3,3,4,4,2)), 2, 3)
  glabs = list("2"='G1', "3"='Early-S', "4"='Late-S/G2')
  for (i in 1:ncol(gind)) {
    g1 = gind[1, i]
    g2 = gind[2, i]

    plot(g_mu[[g1]][valid_doms], g_mu[[g2]][valid_doms], pch=19, cex=cex, xlab=paste(glabs[g1], '(mean)'), ylab=paste(glabs[g2], '(mean)'), main="", xlim=mean_lim, ylim=mean_lim, col=addalpha(pch_col, alpha), asp=1)
    grid(col='black', lwd=0.5*fig_factor)
    abline(a=0, b=1, col='red', lwd=2)

    plot(g_sd[[g1]][valid_doms], g_sd[[g2]][valid_doms], pch=19, cex=cex, xlab=paste(glabs[g1], '(sd)'), ylab=paste(glabs[g2], '(sd)'), main="", xlim=sd_lim, ylim=sd_lim, col=addalpha(pch_col, alpha), asp=1)
    grid(col='black', lwd=0.5*fig_factor)
    abline(a=0, b=1, col='red', lwd=2)
  }
  dev.off()
  
  

  mu_b = apply(mu, 1, function(v) { table(cut(v, breaks=freq_breaks)) } )
  mu_b_n = t(mu_b) / colSums(mu_b, na.rm=T)
  #hc = hclust(dist(mu_b_n), method='ward.D2')
  
  #fig_png(sprintf("%s_clus_frac_hm.png", ofn_pref), nrow(mu_b_n),  nrow(mu_b_n) / 4)
  
  #image(pmin(mu_b_n[ hc$order, ], frac_tab_zmax), zlim=c(0, frac_tab_zmax), col=colorRampPalette(colspec)(nshades), xaxt='n', yaxt='n')
  #axis(2, at=seq(0, 1, len=length(freq_breaks)), labels=freq_breaks, las=1, tick=F)
  #axis(1, at=seq(0, 1, len=nrow(mu_b_n)), labels=1:nrow(mu_b_n), tick=F)
  #dev.off()

  # plot sorted freq table per group (2-4)
  for (g in as.character(2:4)) {
    mu_freq[[g]] = mu_freq[[g]][valid_doms, ]
  }

  if (order_doms_by == "mean_all") {
    #dord = valid_doms[ order(mu_b_n[valid_doms, ] %*% (1:ncol(mu_b_n)), decreasing=T)]
    dord = valid_doms[ order(rowMeans(mu[valid_doms, ], na.rm=T), decreasing=T)]
  }
  else if (order_doms_by == "mean_g4") {
    #dord = valid_doms[ order(mu_freq[['4']] %*% (1:ncol(mu_freq[['4']])), decreasing=T) ]
    dord = valid_doms[ order(g_mu[['4']], decreasing=T) ]
  }
  else if (order_doms_by == "trans_A_score") {
    rownames(d) = paste(d$chrom, d$start, sep="_")
    dord = valid_doms[ order(d[valid_doms, 'trans_A_score'], decreasing=T) ]
  }
  else if (order_doms_by == "all_var") {
    dord = head(valid_doms[ order(apply(mu[valid_doms, ], 1, sd, na.rm=T), decreasing=T) ], 50)
  }
  else if (order_doms_by == "kmeans") {
    ag = do.call("cbind", mu_freq[c('2', '3', '4')])
    km = TGLKMeans_wrapper(ag, sprintf("%s/%s_kmeans", sch_rdata_dir, basename(ofn_pref)), n_clusters)
    clus_m = tapply(g_mu[['4']], km$cluster, mean, na.rm=T)

    km_size = km$size[order(clus_m, decreasing=T)]
    clus_lines = cumsum(km_size) / sum(km_size)
    dord = valid_doms[order(clus_m[km$cluster] + 1e-6 * g_mu[['4']], decreasing=T)]
  }
  

  for (g in as.character(2:4)) {
    b = mu_freq[[g]]
    
    fig_png(sprintf("%s_g%s_b%.2f_by_%s_hm.png", ofn_pref, g, freq_bin, order_doms_by), nrow(b)/4,  nrow(b))
    image(pmin(t(b[ dord, ]), frac_tab_zmax), zlim=c(0, frac_tab_zmax), col=colorRampPalette(freq_colspec)(nshades), xaxt='n', yaxt='n')
    axis(1, at=seq(0, 1, len=length(freq_breaks)), labels=freq_breaks, las=2, tick=F)

    if (order_doms_by == "kmeans") {
      abline(h=clus_lines)
    }
    #axis(2, at=seq(0, 1, len=nrow(b)), labels=1:nrow(b), tick=F)
    dev.off()
  }

  plot_legend(sprintf("%s_hm_leg.png", ofn_pref), c(0, frac_tab_zmax), colorRampPalette(freq_colspec)(nshades), "mean A")

  d_ord = d[dord, grep("color", colnames(d), invert=T)]
  d_ord$ord_by_g4 = 1:nrow(d_ord)
  d_ord$pc = rep(n_pc:1, each=ceiling(nrow(d_ord) / n_pc))[1:nrow(d_ord)]
  
  list(mu_freq=mu_freq, d_ord=d_ord, mean_sd_by_slice=mean_sd_by_slice, mu=mu, valid=valid)
}

plot_domains_weighted_a_score_across_groups_diploid <- function()
{
  rp_dip_as = plot_domains_weighted_a_score_across_groups(freq_lim=c(0.68, 0.84), freq_bin=0.02, frac_tab_zmax=0.3)

}

###############
fig2_index_sort_metrics <- function(idir=sch_table_dir, ifns=paste0("Index_plate_", c(1, 3:8), ".txt"), nm_pref=paste0("scell.nextera.1CDES_p", c(1, 3:8), "_"), columns=c(10, 11), col_names=c('hoechst', 'geminin'), width=320, height=220, pch_col=darkblue_col, cex=0.3*fig_factor, plate_nms=paste0('p', c(1, 3:8)), colspec=rev(brewer.pal(n=11, "RdYlBu")[c(1:4, 8:11)]), only_good=T) 
{
  nms = do.call("c", lapply(nm_pref, grep, x=rownames(sch_decay_metrics), value=T))

  
  print(length(nms))
  df = data.frame(ord=sch_decay_metrics[nms, 'ord'], b=NA)
  rownames(df) = nms
  for (s in col_names) {
    df[, s] = NA
  }

  for (i in seq_along(ifns)) {
    x = read.table(sprintf("%s/%s", idir, ifns[i]), header=T, sep="\t")
    x = filter(x, Population == 'P1')

    df[ grep(nm_pref[i], nms, value=T), 'b'] = i
    
    for (s in seq_along(col_names)) {
      df[ paste0(nm_pref[i], x$Well), col_names[s]] = x[, columns[s]]
    }
  }
  df = df[!is.na(df$ord), ]
  print(nrow(df))
  
  cols = qualitative_colors[df$b]
  cols[df$hoechst > 30000 & df$ord > 100 & df$ord < 1254] = 'black'

  df$bad = cols == 'black'

  if (only_good) {
    df = df[intersect(rownames(df), sch_good_cells), ]
  }

  if (!is.null(pch_col)) {
    cols = pch_col
  }
  
  for(s in col_names) {
    fig_png(sprintf("%s/%s/index_sorting_vs_phased_ord_%s_%s.png", sch_fig_dir, rel_odir, s, ifelse(only_good, "good", "all")), width, height)
    plot(df$ord, df[, s], main=sprintf("cor %.2f", cor(df$ord, df[,s], use='pair')), xlab="phased order", ylab=s, col=cols, pch=19, cex=cex)
    grid(col='black', lwd=0.5*fig_factor)
    abline(v=cumsum(table(sch_decay_metrics$group)))
    dev.off()
  }

  fig_png(sprintf("%s/%s/index_sorting_vs_phased_ord_leg.png", sch_fig_dir, rel_odir), width, height)
  plot(1:10, col=NA, xaxt='n', yaxt='n', xlab="", ylab="", main="")
  legend("topleft", legend=plate_nms, col=qualitative_colors[1:length(plate_nms)], pch=19, bty='n')
  dev.off()


  cols = colorRampPalette(colspec)(nrow(sch_decay_metrics))
  fig_png(sprintf("%s/%s/index_sorting_col_by_phased.png", sch_fig_dir, rel_odir), width, height)
  plot(df[, col_names[1]], df[, col_names[2]], main="", xlab=col_names[1], ylab=col_names[2], col=cols[df$ord], pch=19, cex=cex)
  dev.off()

  df[order(df$ord), 'ord_r'] = 1:nrow(df)
  df$facs_score = scale(df$hoechst) + scale(df$geminin)
  df[order(df$facs_score), 'facs_r'] = 1:nrow(df)
  df$g = sch_decay_metrics[rownames(df), 'group']

  
  fig_png(sprintf("%s/%s/index_sorting_ranks_comparison.png", sch_fig_dir, rel_odir), width, width)

  plot(df$ord_r, df$facs_r, pch=19, col=darkblue_col, cex=cex, asp=1, xlab='phasing rank', ylab="Hoechst + Geminin rank", main=sprintf("g>1 cor=%.2f", cor(df[df$g > 1, 'ord_r'], df[df$g > 1, 'facs_r'], method="spearman")))
  abline(v=cumsum(table(df$g)))
  abline(a=0, b=1, lty=1, lwd=2)


  dev.off()

  df
}

##############
fig_s_pc_pc_contact_enrich <- function(size=400, cis_shuffle_tn=paste0(pool_tn, "_shuffle"), nms=c("g1", "early_s", "mid_s_g2"), nm_cols=mini_qualitative_colors[2:4], pc_ifn=sprintf("%s/domain_ord_by_mean_g4_A_assoc.txt", sch_table_dir), min_cis=2e6, vfunc="weighted.sum", rebuild=F, colspec=diver_colspec, zlim=c(-1.5, 1.5), cex=0.3, alpha=0.5, exp_zlim=NULL)
{
  obs_cis = list()
  obs_trans = list()

  tns = paste(pool_tn, "group", 2:4, nms, sep="_")

  for (i in seq_along(tns)) {
    message("calc obs for ", nms[i])
    obs_cis  [[ nms[i] ]] = calc_pc_pc_total_contacts(tns[i], cis=T, pc_ifn=pc_ifn, min_cis=min_cis, vfunc=vfunc, rebuild=rebuild)
    obs_trans[[ nms[i] ]] = calc_pc_pc_total_contacts(tns[i], cis=F, pc_ifn=pc_ifn, vfunc=vfunc, rebuild=rebuild)
  }

  
  message("calc cis exp from ", cis_shuffle_tn)
  cis_exp = calc_pc_pc_total_contacts(cis_shuffle_tn, cis=T, pc_ifn=pc_ifn, min_cis=min_cis, vfunc=vfunc, rebuild=rebuild)
  cis_exp = cis_exp / sum(cis_exp)

  marg_trans = rowSums(Reduce("+", obs_trans)) 
  marg_trans_n = marg_trans  / sum(marg_trans)
  trans_exp = marg_trans_n %*% t(marg_trans_n)

  marg_cis = rowSums(Reduce("+", obs_cis)) 
  marg_cis_n = marg_cis  / sum(marg_cis)
  cis_marg_exp = marg_cis_n %*% t(marg_cis_n)

  save(obs_cis, obs_trans, cis_exp, trans_exp, file=sprintf("%s/pc_pc_mats.rdata", sch_rdata_dir))
  
  .inner_plot_pc_pc_mat = function(suff, m, zlim, colspec=diver_colspec, vfunc="weighted.sum") {
    fig_png(sprintf("%s/%s/fig_pc_pc_contig_enr_%s_%s.png", sch_fig_dir, rel_odir, vfunc, suff), size, size)
    if (is.null(zlim)) {
      zlim = range(m)
    }
    message(sprintf("%s\t%g\t%g", suff, min(m), max(m)))
    par(mar=c(3,3,3,3))
    image(pmin(pmax(m[, ncol(m):1], zlim[1]), zlim[2]), col=colorRampPalette(colspec)(256), zlim=zlim, main=sprintf("%s %g - %g", suff, zlim[1], zlim[2]), xaxt='n', yaxt='n')
    axis(1, at=seq(0, 1, length=nrow(m)), labels=1:nrow(m), tick=F)
    axis(2, at=seq(1, 0, length=nrow(m)), labels=1:nrow(m), tick=F)
    dev.off()
  }

  message("plotting...")
  .inner_plot_pc_pc_mat("cis_exp", log2(cis_exp), exp_zlim, seq_colspec, vfunc)
  .inner_plot_pc_pc_mat("cis_marg_exp", log2(cis_marg_exp), exp_zlim, seq_colspec, vfunc)

  .inner_plot_pc_pc_mat("trans_exp", log2(trans_exp), exp_zlim, seq_colspec, vfunc)

  ocis_n = NULL
  otrans_n = NULL
  
  for (i in seq_along(nms)) {
    ocis = obs_cis[[ nms[i] ]] / sum(obs_cis[[nms[i] ]])
    ocis_n = rbind(ocis_n, as.vector(ocis))
    
    otrans = obs_trans[[ nms[i] ]] / sum(obs_trans[[ nms[i] ]])
    otrans_n = rbind(otrans_n, as.vector(otrans))
    
    .inner_plot_pc_pc_mat(paste0("cis_", nms[i], "_obs"), log2(ocis), exp_zlim, vfunc=vfunc)
    .inner_plot_pc_pc_mat(paste0("cis_", nms[i], "_norm"), log2(ocis / cis_exp), zlim, vfunc=vfunc)
    .inner_plot_pc_pc_mat(paste0("cis_", nms[i], "_marg_norm"), log2(ocis / cis_marg_exp), zlim, vfunc=vfunc)
    .inner_plot_pc_pc_mat(paste0("trans_", nms[i], "_obs"), log2(otrans), exp_zlim, vfunc=vfunc)
    .inner_plot_pc_pc_mat(paste0("trans_", nms[i], "_norm"), log2(otrans / trans_exp), zlim, vfunc=vfunc)
  }

  e_cis = as.vector(cis_exp)
  cis_lim = range(rbind(ocis_n, e_cis))
  
  fig_png(sprintf("%s/%s/test_pc_contig_obs_vs_exp_cis_%s.png", sch_fig_dir, rel_odir, vfunc), size, size)
  plot(e_cis, ocis_n[1,], col=NA, ylim=range(ocis_n), main='cis')
  apply(cbind(nm_cols, ocis_n), 1, function(x) { points(e_cis, x[-1], pch=19, cex=cex, col=addalpha(x[1], alpha)) })
  grid(col='black', lwd=0.5*fig_factor)
  legend("topleft", legend=paste(nms, apply(ocis_n, 1, function(x) { sprintf("%.2f", cor(x, e_cis)) })), col=nm_cols, pch=19)
  dev.off()

  e_trans = as.vector(trans_exp)
  trans_lim = range(rbind(otrans_n, e_trans))
  fig_png(sprintf("%s/%s/test_pc_contig_obs_vs_exp_trans_%s.png", sch_fig_dir, rel_odir, vfunc), size, size)
  plot(e_trans, otrans_n[1,], col=NA, ylim=range(otrans_n), main='trans')
  apply(cbind(nm_cols, otrans_n), 1, function(x) { points(e_trans, x[-1], pch=19, cex=cex, col=addalpha(x[1], alpha)) })
  grid(col='black', lwd=0.5*fig_factor)
  legend("topleft", legend=paste(nms, apply(otrans_n, 1, function(x) { sprintf("%.2f", cor(x, e_trans)) })), col=nm_cols, pch=19)
  dev.off()
}

#################
fig_s_domains_epi_stats <- function(size=400, pc_ifn=sprintf("%s/domain_ord_by_mean_g4_A_assoc.txt", sch_table_dir), cex=0.3*fig_factor, trim_q=0.002, trimmed_col=NA)
{
  pc = read.table(pc_ifn, header=T)
  pc = pc[order(pc$ord_by_g4, decreasing=T), ]
  
  exprs = list(RNA=sprintf("log2((%s + %s)/2)", sch_rna_tn_rep1, sch_rna_tn_rep2),
    H3K4me1="(Encode.esb4.h3k4me1.rep1 + Encode.esb4.h3k4me1.rep2)/2",
    H3K4me3="(Encode.esb4.h3k4me3.rep1 + Encode.esb4.h3k4me3.rep2)/2",
    laminB1="laminB1.ESC")

  for (s in names(exprs)) {
    x = gextract(exprs[[s]], intervals=pc, iterator=pc, colnames=s)

    ylim = quantile(x[, s], c(trim_q, 1 - trim_q), na.rm=T)
    
    outliers = x[, s] < ylim[1] | x[, s] > ylim[2]
    
    fig_png(sprintf("%s/%s/fig4_%s_by_A_assoc.png", sch_fig_dir, rel_odir, s), size, size)
    
    plot(x$intervalID, pmin(pmax(x[, s], ylim[1]), ylim[2]), pch=19, cex=cex, col=ifelse(outliers, trimmed_col, darkblue_col), xlab='ord domains', ylab=s, main="")
    grid(col='black', lwd=0.5*fig_factor)
    dev.off()
  }
    
    
}


# Plot genome-wide idiogram of domains colored by given colors
sch_plot_domain_score_idiogram <- function(domains, dom_grouping, fn, colors, legend_nms, chrs=paste0('chr', sch_chroms), chr_space=0.1, title="")
{
    chroms = gintervals.all()
    rownames(chroms) = chroms$chrom

    chroms = chroms[chrs, ]
    chroms$ord = 1:nrow(chroms)

    height = 36 * nrow(chroms)
    width  = height * (4/3)

    fig_png(fn, width=width, height=height)

    plot.new()
    plot.window(xlim=c(0, max(chroms$end)), ylim=c(0,nrow(chroms)))

    domains$ord = chroms[as.vector(domains$chrom), "ord"]

    y0 = nrow(chroms)
    rect(xleft=domains$start, xright=domains$end, ybottom=y0 - (chr_space + domains$ord - 1), ytop=y0 - (domains$ord - chr_space), col=colors[dom_grouping], border=NA)
    # rect(xleft=domains$start, xright=domains$end, ybottom=(chr_space + domains$ord - 1), ytop=(domains$ord - chr_space), col=colors[dom_grouping], border=NA)
    rect(xleft=0, xright=chroms$end, ybottom=rev(seq(chr_space, by=1, length=nrow(chroms))), ytop=rev(seq(1-chr_space, by=1, length=nrow(chroms))), col=NA, border="black")

    text(x=0, y=rev(seq(0.5, by=1, length=nrow(chroms))), labels=gsub(chrs, pattern="chr", replace=""), pos=2, cex=1.5)
    title(main=title, cex.main=2)
    if (!is.null(legend_nms)) {
        legend("bottomright", legend_nms, col=colors, pch=15, bty="n", inset=c(0,0), cex=1.5)
    }
    dev.off()
}
