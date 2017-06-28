###############################################
#
# Clustering/Ordering of cells
#
###############################################


#===================
plot_clus_decays_and_stats <- function(mean_cmap_rel_dir="cmaps/far_n/", clus=NULL, width=120, height=120, other_col='black', clus_col='blue', pcex=0.1, only_decays=F, only_scatters=F)
{
  if (is.null(clus)) {
    clus = sch_decay_clust
  }
  k = max(clus)

  x = sch_glob_decay_res[sch_good_cells, -1]
  mv = sch_decay_metrics[sch_good_cells, ]
  mv$clus = clus$x

  x = x / rowSums(x)
  x = x[, 1:(ncol(x)-1)]

  col_dists = as.numeric(gsub("X", "", colnames(x)))

  dc = list()
  ylims = 0
  for (i in 1:k) {
    dc[[i]] = colMeans(x[ rownames(x)[clus == i], ])
    ylims = c(ylims, range(dc[[i]]))
  }

  for (i in 1:k) {
    if (only_decays) {
      png(sprintf("%s/%s/c%d_line_decay.png", sch_fig_dir, mean_cmap_rel_dir, i), width=width, height=height*1.5, pointsize=pointsize, res=fres)
      par(mar=c(4,3,0,0))
    }
    else if (only_scatters) {
      png(sprintf("%s/%s/c%d_scatters.png", sch_fig_dir, mean_cmap_rel_dir, i), width=width, height=height*2.5, pointsize=pointsize, res=fres)
      layout(matrix(1:2, 2, 1))
      par(mar=c(2,4,0,0))
    }
    else {
      png(sprintf("%s/%s/c%d_line_decay.png", sch_fig_dir, mean_cmap_rel_dir, i), width=width, height=height*3, pointsize=pointsize, res=fres)
      layout(matrix(1:3, 3, 1))
      par(mar=c(5,3,0,0))
    }

    if (!only_scatters) {
      plot(log2(col_dists), dc[[i]], type='l', ylab="", xlab=sprintf("C%d (%d)", i, sum(clus == i)), ylim=range(ylims), main="", lwd=2, col=clus_col)
      grid(col='black', lwd=0.5)
      for (j in 1:k) {
        lines(log2(col_dists), dc[[j]], lwd=1, col=other_col)
      }
      lines(log2(col_dists), dc[[i]], lwd=2, col=clus_col)


      #inds = seq(1, ncol(x), by=3)
#    axis(1, at=inds, labels=sapply(col_dists[inds], n2str), cex=0.5)
      #text(x=log2(col_dists[1]), y=0.9*max(ylims), labels=sprintf("C%d (%d)", i, sum(clus == i)), pos=4)
    }

    if (!only_decays) {
      #plot(mv$high_near_f, mv$far_mu, main="", xlab="", ylab="", pch=19, cex=pcex, col=other_col)
      plot(mv$f_near_band, mv$far_mu, main="", xlab="", xaxt='n', ylab="", pch=19, cex=pcex, col=other_col)
      ccells = rownames(x)[clus == i]
      points(mv[ccells, 'f_near_band'], mv[ccells, 'far_mu'], pch=19, cex=pcex, col=clus_col)
      grid(col='black', lwd=0.5)

      mv$early_f = mv$early_f / min(mv$early_f, na.rm=T)
      plot(mv$f_near_band, mv$early_f, main="", ylab=sprintf("C%d (%d)", i, sum(clus == i)), xlab="", pch=19, cex=pcex, col=other_col)
      points(mv[ccells, 'f_near_band'], mv[ccells, 'early_f'], pch=19, cex=pcex, col=clus_col)
      grid(col='black', lwd=0.5)
    }

    dev.off()


  }
  print(sort(tapply(mv$f_near_band, mv$clus, mean)))
}

#############
# clus_order=c(1,3,10,9,13,2,11,8,12,5,6,4,7))
plot_decay_mat_sorted_by_x <- function(bin="X11863283.2030314", zlim=c(0, 0.035), ec_zlim=c(0.45, 0.75), ab_zlim=c(0.15, 0.95), order_by="comp_cluster_peak_ord", k=30, clus_order="hclust", intra_clus_ord_func=mean, clus_ord_decr=F, which_max_smooth_width=5, domain_peak_dists=2**c(11.25, 21.75), far_cis_dists=2**c(23, 27.25), chroms=NULL, name="", remove_mitotic=F, mitotic_dists=2**c(21, 22), mitotic_mean_threshold=0.016, use_cells=NULL, hclust_clus_by="euclid", remove_trans=F, cells_names_ord=NULL, matrix_type="all", add_dist_grid=T, add_ab=T, add_nodig_trans=F, remove_unear_dist=NULL, vlines=NULL, width=NULL, height=NULL, ext_ofn=NULL, bins=2**c(21, 23.5), plot_lines=F, only_hm=F, colspec=c('blue', 'white', 'red'))
{
  if (is.null(chroms)) {
    x = switch(matrix_type,
      all=sch_glob_decay_res,
      early=sch_tor_glob_decay[['e']],
      mid_early=sch_tor_glob_decay[['me']],
      mid_late=sch_tor_glob_decay[['ml']],
      late=sch_tor_glob_decay[['l']])
    x = x[sch_good_cells, -1]
  }
  else {
    stopifnot(matrix_type == "all")
    x = sch_chrom_decay_res[sch_good_cells,]
    cols = c()
    for (chr in chroms) {
      cols = c(cols, grep(paste("chr", chr, "_", sep=""), colnames(x), value=T))
    }
    x = x[,cols]
    x = t(apply(x, 1, tapply, gsub("chr.*_", "X", cols, perl=T), sum))
    colnames(x) = gsub("Xtrans", "trans", colnames(x))
    x = x[, colnames(sch_glob_decay_res)[-1]]
  }

  if (!is.null(use_cells)) {
    x = x[ intersect(rownames(x), use_cells), ]
  }

  if (remove_trans) {
    x = x[, -ncol(x)]
  }

  if (!is.null(remove_unear_dist)) {
    next_unear_bin = min(which(as.numeric(gsub("X", "", colnames(x))) > remove_unear_dist))
    message("starting from ", next_unear_bin)
    x = x[, next_unear_bin:ncol(x)]
  }
  x = x / rowSums(x)

  if (remove_mitotic) {
    stopifnot(is.null(remove_unear_dist))
    mitos_bins = which (sch_cis_breaks >= mitotic_dists[1] & sch_cis_breaks <= mitotic_dists[2]) - 2
    message ("removing ",  sum(apply(x[, mitos_bins], 1, mean) >= mitotic_mean_threshold), " mitotic cells")
    x = x[apply(x[, mitos_bins], 1, mean) < mitotic_mean_threshold, ]
    print(dim(x))
  }

  ofn_pref = sprintf("%s/ordering/decay_clust_sort_%s_%s%s%s%s%s", sch_fig_dir, matrix_type, name, switch(order_by, which_max="which_max", far_near_ratio="far_near_ratio", early_cov="early_cov", peaks_ratio="peaks_ratio", bin=paste(order_by, bin, sep="_"), bins=paste(sapply(bins, n2str, 2), collapse="_"), comp_cluster=paste("comp_cluster_k", k, sep=""), comp_cluster_peak_ord=paste("comp_cluster_peaks_ord_k", k, sep=""),  comp_cluster_farcis_ord=paste("comp_cluster_farcis_ord_k", k, sep=""), comp_cluster_hclust=paste(order_by, k, "_", hclust_clus_by, sep=""), hclust="hclust", cis_aa_enrich="cis_aa_enrich", cis_bb_enrich="cis_bb_enrich", given="given"), ifelse(remove_mitotic, "_noM", ""), ifelse(remove_trans, "_noTrans", ""), ifelse(is.null(remove_unear_dist), "", sprintf("_use_gt%s", n2str(remove_unear_dist))))

  height = ifelse(is.null(height), 800, height)
  width = ifelse(is.null(width), max(1500, nrow(x)+100), width)

  fig_png(ifelse(is.null(ext_ofn), paste(ofn_pref, ".png", sep=""), ext_ofn), height=height, width=width)

  if (only_hm) {
    lmat = matrix(1, 1, 1)
    lheights = 1
  }
  else {
    if (add_ab) {
      if (add_nodig_trans) {
        lmat = matrix(6:1, 6, 1)
        lheights = c(1,1,3,1,1,14)
      }
      else {
        lmat = matrix(4:1, 4, 1)
        lheights = c(3,1,1,14)
      }
    }
    else {
      if (add_nodig_trans) {
        lmat = matrix(5:1, 5, 1)
        lheights = c(1,1,3,1,14)
      }
      else {
        lmat = matrix(3:1, 3, 1)
        lheights = c(3,1,14)
      }

    }
  }
  if (plot_lines) {
    lheights[length(lheights)-1] = lheights[length(lheights)-1] * 4
    lheights[length(lheights)] = lheights[length(lheights)] * 2
  }

  print(lheights)
  layout(lmat, heights=lheights)

  par(mar=c(1,5,1,1))

  if (add_ab) {
    ab = calc_ab_enrich(sch_ab_counts[rownames(x),  grep("cis", colnames(sch_ab_counts))]/2, "cis")
  }

  xdists = as.numeric(gsub("X", "", colnames(x)))

  if (is.null(cells_names_ord)) {
    if (length(grep("comp_cluster", order_by)) > 0) {
      clus_order = if (order_by == "comp_cluster") clus_order else NULL
      order_by_str = switch(order_by,
        comp_cluster=NULL,
        comp_cluster_peaks_ord="peaks_ratio",
        comp_cluster_farcis_ord="farcis_frac",
        comp_cluster_hclust="hclust")

      km = cluster_decay_profiles(k=k, clus_order=clus_order, intra_clus_ord_func=intra_clus_ord_func, clus_ord_decr=clus_ord_decr, order_by_str=order_by_str, remove_mitotic=remove_mitotic, domain_peak_dists=domain_peak_dists, far_cis_dists=far_cis_dists, hclust_clus_by=hclust_clus_by, remove_trans=remove_trans, mitotic_mean_threshold=mitotic_mean_threshold)
      cells_ord = km$order
    }
    else if (order_by == "which_max") {
                                        #ord = t(apply(x[,-ncol(x)], 1, order, decreasing=T))
                                        #cells_ord = order(ord[,1], ord[,2], ord[,3], ord[,4], ord[,5], ord[,6], ord[,7], ord[,8], ord[,9], ord[,10], decreasing=T)
      cells_ord = order(apply(x[,-ncol(x)], 1, function(v) { which.max(rollapply(v, width=which_max_smooth_width, by=1, FUN=mean, partial=T)) }))
    }
    else if (order_by == "bin") {
      cells_ord = order(x[,bin], decreasing=T)
    }
    else if (order_by == "bins") {
      cells_ord = order(rowSums(x[, xdists > bins[1] & xdists <= bins[2]]), decreasing=T)
    }
    else if (order_by == "far_near_ratio") {
      cells_ord = order(rowSums(x[,3:49])/rowSums(x[50:69]))
    }
    else if (order_by == "early_cov") {
      cells_ord = order(rowSums(sch_tor_marg_n[rownames(x), c(3,4)]))
    }
    else if (order_by == "peaks_ratio") {
      cells_ord = order(apply(x, 1, get_peaks_ratio, p1=domain_peak_dists, p2=far_cis_dists))
    }
    else if (order_by == "cis_aa_enrich") {
      cells_ord = order(ab$enAA)
    }
    else if (order_by == "cis_bb_enrich") {
      cells_ord = order(ab$enBB)
    }
    else if (order_by == "hclust") {
      hc = hclust(dist(x[, 16:65]), method="ward.D2")
      dd = as.dendrogram(hc)
      reord.dd = reorder(dd, apply(x, 1, get_peaks_ratio, p1=domain_peak_dists, p2=far_cis_dists), mean)
      reord.hc = as.hclust(reord.dd)
      cells_ord = reord.hc$order
    }
    else {
      message(sprintf("unknown order_by parameter: %s", order_by))
      return(NA)
    }
  }
  else {
    x = x[cells_names_ord, ]
    cells_ord = 1:nrow(x)
  }

  image(pmin(pmax(as.matrix(x[cells_ord,]), zlim[1]), zlim[2]), col=colorRampPalette(colspec)(100), zlim=zlim, xaxt='n', yaxt='n', xlab="", ylab="", main="")
  sbins = sapply(as.numeric(gsub("X", "", colnames(x))), n2str, 2)

  ax_at = seq(0, 1, len=ncol(x))
  ax_c = seq(1, ncol(x), by=4)
  if (!remove_trans) {
    sbins[length(sbins)] = 'trans'
  }
  axis(2, ax_at[ax_c], sbins[ax_c], las=2, cex.axis=1.4)

  if (order_by == "comp_cluster" || order_by == "comp_cluster_peak_ord") {
    abline(v=cumsum(km$size)/sum(km$size), lwd=2)
  }

  if (add_dist_grid) {
    grid_dists = sort(unique(as.vector(c(sch_far_dists, sch_high_near_dists))))
    grid_cols = which(is.element(colnames(x), paste0("X", grid_dists)))
    abline(h=grid_cols / ncol(x), col='black')
  }

  if (!is.null(vlines)) {
    abline(v=vlines, lwd=2)
  }

  if (order_by == "bins") {
    abline(h=range(which(xdists > bins[1] & xdists <= bins[2]))/length(xdists))
  }

  if (!only_hm) {
    yy = rowSums(sch_tor_marg_n[rownames(x)[cells_ord],c(3,4)])
    cov_col = colorRampPalette(c("green", "black", "red"))(200)
    if (plot_lines) {
      z = (yy - min(yy)) / (max(yy) - min(yy))
      image(matrix(0, length(z), length(z)), col='white', axes=F, main="Early cov")
                                        #points(z, seq(0, 1, length=length(z)), pch=19, col=cov_col[1 + ceiling(z*(length(cov_col)-1))])
      points(z, seq(0, 1, length=length(z)), pch=19)
      axis(1)
      box()
      grid(col='darkgray')
    }
    else {
      image(pmin(pmax(as.matrix(yy), ec_zlim[1]), ec_zlim[2]), col=colorRampPalette(c("green", "black", "red"))(200), zlim=ec_zlim, xaxt='n', yaxt='n', xlab="", ylab="ecov", main="")
    }
    message(sprintf("early cov: trimming %f - %f to %f - %f", min(yy), max(yy), ec_zlim[1], ec_zlim[2]))


                                        #  image(pmin(pmax(as.matrix(ab[rownames(x)[cells_ord], 'enAA']), 0.15), 0.55), col=colorRampPalette(c("blue", "black", "yellow"))(200), zlim=c(0.15, 0.55), xaxt='n', yaxt='n', xlab="", ylab="", main="")
                                        #  image(pmin(pmax(as.matrix(ab[rownames(x)[cells_ord], 'enBB']), 0.1), 0.8), col=colorRampPalette(c("blue", "black", "yellow"))(200), zlim=c(0.1, 0.8), xaxt='n', yaxt='n', xlab="", ylab="", main="")
    if (add_ab) {
      yy = ab[rownames(x)[cells_ord], 'deplAB']
      image(pmin(pmax(as.matrix(yy), ab_zlim[1]), ab_zlim[2]), col=colorRampPalette(c("blue", "black", "yellow"))(200), zlim=ab_zlim, xaxt='n', yaxt='n', xlab="", ylab="ab", main="")
      message(sprintf("ab: trimming %f - %f to %f - %f", min(yy), max(yy), ab_zlim[1], ab_zlim[2]))
    }

                                        #  image(pmin(pmax(as.matrix(t(rowMeans(sch_fend_dup[rownames(x)[cells_ord],]))),0),0.02),         col=colorRampPalette(c("white", "brown"))(100), zlim=c(0, 0.02), xaxt='n', yaxt='n', xlab="", ylab="", main="")

    cond = sch_batch[rownames(x), 'cond']
    batch_col = c("white", unlist(sch_cond_colors))

    cb = matrix(0, length(cond), length(sch_cond_colors))
    colnames(cb) = names(sch_cond_colors)
    rownames(cb) = rownames(x)
    cond_idx = seq_along(sch_cond_colors)
    names(cond_idx) = names(sch_cond_colors)
    cb[ cbind(rownames(cb), cond)] = cond_idx[cond]

                                        #print(cb[cells_ord,])
    image(as.matrix(cb[cells_ord,]), col=batch_col, xaxt='n', yaxt='n', xlab="", ylab="", main="")
    axis(2, seq(0, 1, len=length(cond_idx)), labels=names(cond_idx), las=2)


    if (add_nodig_trans) {
      dn = sch_glob_decay_res[sch_good_cells,]
      dn = dn / rowSums(dn)
      image(as.matrix(dn[rownames(x)[cells_ord], 1]), col=colorRampPalette(c("white", "black"))(100), xaxt='n', yaxt='n', xlab="", ylab="", main="")

      image(as.matrix(dn[rownames(x)[cells_ord], ncol(sch_glob_decay_res)]), col=colorRampPalette(c("white", "black"))(100), xaxt='n', yaxt='n', xlab="", ylab="", main="")
    }

    ## ord_clus = rownames(x)[cells_ord]
    ## cm = matrix(0, nrow(x), sch_decay_k)
    ## cm[ cbind(seq_along(ord_clus), sch_decay_clust[ord_clus,])] = 1
    ## image(as.matrix(t(cm)), col=c("white", "black"), xaxt='n', yaxt='n', xlab="", ylab="", main="")
    ## axis(1, at=seq(0, 1, length=sch_decay_k), labels=1:sch_decay_k)
  }
  dev.off()

  list(ofn=ofn_pref, cells=rownames(x)[cells_ord])
}

##
#get_peaks_ratio=function(v, p1=2**c(15.25, 20.5), p2=2**c(23.5, 26.75)) {
get_peaks_ratio=function(v, p1=2**c(11.25, 21.75), p2=2**c(23, 27.25)) {
  p1_bins = paste("X", sch_cis_breaks[which (sch_cis_breaks >= p1[1] & sch_cis_breaks <= p1[2])], sep="")
  p2_bins = paste("X", sch_cis_breaks[which (sch_cis_breaks >= p2[1] & sch_cis_breaks <= p2[2])], sep="")
  sum(v[p1_bins], na.rm=T)/sum(v[p2_bins], na.rm=T)
}




#############
# 1. First cluster - mitotic outliers (cluster #1)
# 2. Intra cluster order by which.max bin
# 3. Allow to order clusters by a given permutation (1:k, 1 is mitotic)
cluster_decay_profiles <- function(mitotic_dists=2**c(21, 22), mitotic_mean_threshold=0.016, k=12, cluster_dists=2**c(11.25, 27.25), clus_order=NULL, remove_no_digest=T, remove_trans=F, intra_clus_ord_func=which.max, clus_ord_decr=F, order_by_str=NULL, remove_mitotic=F, domain_peak_dists=2**c(15.25, 20.5), far_cis_dists=2**c(23.5, 26.75), hclust_clus_by="pearson")
{
  x = sch_glob_decay_res[sch_good_cells, ]
  if (remove_no_digest) {
    x[,1] = 0
  }
  if (remove_trans) {
    x[,ncol(x)] = 0
  }

  x = x / rowSums(x)

  clus = rep(1, nrow(x))
  names(clus) = sch_good_cells

  # remove mitotic cells
  mitos_bins = which (sch_cis_breaks >= mitotic_dists[1] & sch_cis_breaks <= mitotic_dists[2]) - 1
  y = x[apply(x[, mitos_bins], 1, mean) < mitotic_mean_threshold, ]

  # cluster
  cluster_bins = which (sch_cis_breaks >= cluster_dists[1] & sch_cis_breaks <= cluster_dists[2]) - 1
  km = TGLKMeans_wrapper(y[, cluster_bins], "/tmp/decay_profile_clus.tab", k)
  clus[ names(km$cluster) ] = km$cluster + 1

  if (remove_mitotic) {
    clus = clus[ clus > 1] - 1
    x = y
  }

  order_by_vec = apply(x[, cluster_bins], 1, intra_clus_ord_func)

  if (!is.null(clus_order)) {
    clus2 = clus_order[clus]
    names(clus2) = names(clus)
    clus = clus2
  }
  else if (!is.null(order_by_str)) {
    if (order_by_str == "peaks_ratio") {
      order_by_vec = apply(x, 1, get_peaks_ratio, p1=domain_peak_dists, p2=far_cis_dists)
      clus_order = order(tapply(order_by_vec, clus, intra_clus_ord_func))
    }
    else if (order_by_str == "farcis_frac") {
      p2 = far_cis_dists
      p2_bins = paste("X", sch_cis_breaks[which (sch_cis_breaks >= p2[1] & sch_cis_breaks <= p2[2])], sep="")
      order_by_vec = rowSums(x[, p2_bins])
      clus_order = order(tapply(order_by_vec, clus, intra_clus_ord_func))

    }
    else if (order_by_str == "hclust") {
      centers = km$centers
      if (!remove_mitotic) {
        centers = rbind(colMeans(x[clus == 1, cluster_bins]), centers)
      }
      colnames(centers) = colnames(x)[cluster_bins]

      order_by_vec = apply(x, 1, get_peaks_ratio, p1=domain_peak_dists, p2=far_cis_dists)

      cdist = switch(hclust_clus_by,
        euclid=dist(centers),
        pearson=as.dist(1 - cor(centers)))

      hc = hclust(cdist)
      dd = as.dendrogram(hc)

      reord.dd = reorder(dd, apply(centers, 1, get_peaks_ratio, p1=domain_peak_dists, p2=far_cis_dists), mean)
      reord.hc = as.hclust(reord.dd)
      clus_order = reord.hc$order
    }
    clus2 = order(clus_order)[clus]

    names(clus2) = names(clus)
    clus = clus2
  }

  km$order = order(clus + ifelse(clus_ord_decr, -1, 1) * 1e-5 * order_by_vec)

  km$size = table(clus)
  km
}



###################
#
# V2 stuff
#
###################
#                  haploid (no rem unear)   hyb_es (unear=2**14.75)
# dom_cutoff       0.45 (0.7 when removing unear)
# far_tight_cutoff 0.06 (0.053)
# near_f_cuttoff   0.85 (0.8 when removing unear)
# sch_glob_min_wm_dist_bin      bin number 40 (8 when removing unear)
#
# Hyb (before major X1-4 batches): early_cp=c("green", "black", "red"), k=30, dom_cutoff=0.7, far_tight_cutoff=0.053, near_f_cuttoff=0.8, mitotic_cutoff=0.0105, mitotic_bin= 2**20.875, add_ab=F, add_nodig_trans=F, plot_extra=F)
# prev mitotic_cutoff = 0.01175
#
# hap:
# mitotic > 0.013
# high_near_f cutoofs: 0.66 0.82
group_cells_by_metrics <- function(early_cp=c("green", "black", "red"), k=30, dom_cutoff=0.69, far_tight_cutoff=0.053, near_f_cuttoff=0.8, add_ab=F, add_nodig_trans=F, plot_basic=T, plot_extra=F, criteria="newer", cells=rownames(sch_decay_metrics), chroms=NULL, use_partitions=NULL, n_partitions=0)
{
  stopifnot(criteria %in% c("orig", "new", "newer"))
  use_cells = intersect(sort(cells), rownames(sch_batch)[is.element(sch_batch$batch_id, sch_use_batches)])

  
  if (is.null(chroms)) {
    if (is.null(use_partitions)) {
      dec = sch_glob_decay_res[use_cells, ]
      chrom_to_glob = F
    }
    else {
      x = Reduce("+", lapply(use_partitions, function(i) { read.table(sprintf("%s/chrom_decay_res_step%s_p%d_of%d.txt", sch_table_dir, as.character(sch_decay_step), i, n_partitions), header=T, sep="\t") } ))
      x = x[use_cells, ]
      chrom_to_glob=T
    }
  }
  else {
    chroms = gsub("chr", "", chroms)
    x = as.matrix(sch_chrom_decay_res[use_cells, grep(paste(paste0("chr", chroms, "_"), collapse="|"), colnames(sch_chrom_decay_res), perl=T)])
    chrom_to_glob = T
  }

  if (chrom_to_glob) {
    xc = t(apply(x, 1, function(v) { tapply(X=v, INDEX=gsub("chr.*_", "", colnames(x), perl=T), FUN=sum) }))
    colnames(xc) = gsub("Xtrans", "trans", paste0("X", colnames(xc)))
    dec = xc[, colnames(sch_glob_decay_res)]
  }
  #print(mean(rowSums(dec)))
  
  dec_dists = c(sch_cis_breaks[-1], -1)

  last_to_remove = 1
  if (!is.null(sch_remove_near_dist)) {
     last_to_remove = max(which(dec_dists > 0 & dec_dists <= sch_remove_near_dist))
     #message("removing 1 to ", last_to_remove)
  }

  dind = -c(1:last_to_remove, ncol(dec))
  dec_n = dec[, dind]
  dec_n = dec_n / rowSums(dec_n)
  dec_n_dists = dec_dists[dind]
  

  all_mu = as.matrix(dec_n) %*% dec_dists[-c(1:last_to_remove, ncol(dec))]

  f_mit_band    = rowSums(dec_n[, dec_n_dists > sch_mitotic_dists[1] & dec_n_dists <= sch_mitotic_dists[2]])
  f_near   = rowSums(dec_n[, dec_n_dists > sch_near_dists[1] & dec_n_dists <= sch_near_dists[2]])

  near_cols = dec_dists >= sch_near_dists[1] & dec_dists <= sch_near_dists[2]
  near = dec[, near_cols]
  near_dists = dec_dists[near_cols]
  near_n = near / rowSums(near)

  domain_cols = near_dists >= sch_domain_dists[1] & near_dists <= sch_domain_dists[2]
  domain = dec[, domain_cols]
  domain_dists = dec_dists[domain_cols]
  domain_n = domain / rowSums(domain)
  #domain_mu = as.matrix(domain_n) %*% domain_dists

  high_near_cols = dec_dists >= sch_high_near_dists[1] & dec_dists <= sch_high_near_dists[2]
  high_near = dec[, high_near_cols]
  high_near_dists = dec_dists[high_near_cols]
  high_near_n = high_near / rowSums(high_near)
  high_near_mu = as.matrix(high_near_n) %*% high_near_dists

  mitotic_cols = dec_dists >= sch_mitotic_dists[1] & dec_dists <= sch_mitotic_dists[2]
  mitotic = dec[, mitotic_cols]
  mitotic_dists = dec_dists[mitotic_cols]
  mitotic_n = mitotic / rowSums(mitotic)

  far_cols = dec_dists >= sch_far_dists[1] & dec_dists <= sch_far_dists[2]
  far = dec[, far_cols]
  far_dists = dec_dists[far_cols]
  far_n = far / rowSums(far)
  far_mu = as.matrix(far_n) %*% far_dists

  #low_near_cols = (dec_dists >= sch_low_near_dists[1] & dec_dists <= sch_low_near_dists[2])[c(-1, -ncol(dec))]
  domain_cols   = (dec_dists >= sch_domain_dists[1] & dec_dists <= sch_domain_dists[2])[c(-1, -ncol(dec))]

  m = data.frame(
    far_tightness = apply(far_n, 1, max),
    far_wm = apply(far_n, 1, which.max), 
    far_mu=far_mu,
    near_wm = apply(near_n, 1, which.max),
    dom_score = rowSums(near_n[, near_dists >= sch_domain_dists[1] & near_dists <= sch_domain_dists[2]]),
    near_f = rowSums(near) / rowSums(dec[, near_cols | far_cols]),
    high_near_f = rowSums(high_near) / rowSums(dec[, high_near_cols | far_cols]),
    high_near_st = rowSums(high_near) / rowSums(dec),
    high_near_mu=high_near_mu,
    #domain_mu=domain_mu,
    all_mu=all_mu,
    mitotic_bin=dec_n[, paste0("X", sch_mitotic_bin)],
    glob_max=log2(as.numeric(gsub("X", "", colnames(dec_n)[ apply(dec_n, 1, which.max)]))),
    f_mitotic_band=f_mit_band,
    f_near_band=f_near)

  m$valid = m$glob_max > log2(sch_glob_min_wm_dist_bin)
  #print(colnames(dec_n)[glob_min_wm])

  #print(quantile(m$mitotic_bin, seq(0,1, by=0.01)))
  #message(sum(m$mitotic_bin > sch_mitotic_cutoff), " mitotic cells")
  #message(sprintf("valid: %d (%d removed)", sum(m$valid), sum(!m$valid)))

  rownames(m) = rownames(dec)

  m$batch = sch_batch[rownames(dec), 'batch']
  m$cond  = sch_batch[rownames(dec), 'cond']

  m$colors = unlist(sch_cond_colors[m$cond])

  if (is.null(use_partitions)) {
    tor_stats = sch_tor_stat
  }
  else {
    tor_stats = Reduce("+", lapply(use_partitions, function(i) { read.table(sprintf("%s/tor_stat_p%d_of%d.txt", sch_table_dir, i, n_partitions), header=T, sep="\t") } ))    
  }

  if (!is.null(chroms)) {
    tor_stats = tor_stats[, is.element(gsub(".*_.*_chr", "", colnames(tor_stats)), chroms)]
  }

  tor_marg = t(apply(tor_stats, 1, function(v) { tapply(v, gsub(".*_tor", "tor", gsub("_chr.*", "", colnames(tor_stats))), sum) } ))
  
  #print(mean(rowSums(tor_marg)))
  
  tor_marg_n =  tor_marg / rowSums(tor_marg)
  m$early_f = rowSums(tor_marg_n[rownames(m), c('tor3', 'tor4')])
  

  m$ec = colorRampPalette(early_cp)(200)[vals_to_cols(m$early_f, seq(min(m$early_f), max(m$early_f), len=3), 100)]

  mv = m[m$valid,]
  
  if (plot_basic) {
    system(sprintf("mkdir -p %s/ordering", sch_fig_dir))
    png(sprintf("%s/ordering/metrics_bp.png", sch_fig_dir), 1000, 180)
    layout(matrix(1:6, 1, 6))
    par(mar=c(3, 4, 3, 1))

    plot(density(m$far_tightness), lwd=2, xlab="", main="far_tightness")
    grid(col='black')

    hist(m$far_wm, ncol(far_n), col='blue', border=NA, main="far_wm", xlab="")
    grid(col='black')

    plot(density(m$dom_score), lwd=2, xlab="", main="domain score")
    grid(col='black')

    hist(m$near_wm, ncol(near_n), col='blue', border=NA, xlab="", main="near_wm")
    grid(col='black')

    plot(density(m$near_f), lwd=2, xlab="", main="frac_near")
    grid(col='black')

    hist(log2(as.numeric(gsub("X", "", colnames(dec_n)))[apply(dec_n, 1, which.max)]), ncol(dec_n), col='blue', border=NA, main="glob wm", xlab="")
    abline(v=log2(sch_glob_min_wm_dist_bin), col='red')
    grid(col='black')

    dev.off()

    # scatters
    .scatter_metric_pairs = function(x, y, lx, ly, cex=0.6, cex.lab=1.5) {
      plot(x, y, pch=19, cex=cex, col=m$colors, main="", xlab=lx, ylab=ly, cex.lab=cex.lab)
      grid()
      plot(x, y, pch=19, cex=cex, col=m$ec, main="", xlab=lx, ylab=ly, cex.lab=cex.lab)
      grid()
    }

    sn = c("dom_score", "high_near_f", "near_wm", "far_tightness", "far_wm", "far_mu", "high_near_mu")

    png(sprintf("%s/ordering/metrics_comp_scatters.png", sch_fig_dir), 1200, 500)
    layout(matrix(1:10, 2, 5))
    par(mar=c(5, 4, 1, 1))
    .scatter_metric_pairs(m[,sn[1]], m[,sn[4]], sn[1], sn[4])
    .scatter_metric_pairs(m[,sn[1]], m[,sn[2]], sn[1], sn[2])
    .scatter_metric_pairs(m[,sn[1]], m[,sn[5]] + 20 * m[,sn[4]], sn[1], sprintf("%s + 20 * %s", sn[5], sn[4]))
    .scatter_metric_pairs(m[,sn[2]], m[,sn[4]], sn[2], sn[4])
    .scatter_metric_pairs(m[,sn[3]] + 20 * m[,sn[1]], m[,sn[5]] + 20 * m[,sn[4]], sprintf("%s + 20 * %s", sn[3], sn[1]), sprintf("%s + 20 * %s", sn[5], sn[4]))

    dev.off()

    # clustering
    tfn = sprintf("/tmp/stam%s", Sys.getpid())
    kf = TGLKMeans_wrapper(far_n, tfn, k)
    kn = TGLKMeans_wrapper(near_n, tfn, k)

    m$kfc = kf$cluster
    m$knc = kn$cluster

    fc = colorRampPalette(early_cp)(200)[vals_to_cols(tapply(m$early_f, m$kfc, mean), seq(min(m$early_f), max(m$early_f), len=3), 100)]
    nc = colorRampPalette(early_cp)(200)[vals_to_cols(tapply(m$early_f, m$knc, mean), seq(min(m$early_f), max(m$early_f), len=3), 100)]

    png(sprintf("%s/ordering/decay_clusters_k%d.png", sch_fig_dir, k), 900, 500)
    layout(matrix(1:2, 1, 2))
    par(mar=c(5,3,3,1))

    plot(colMeans(near_n), type='l', lwd=4, ylim=range(rollapply(kn$centers, width=5, mean, partial=T)), main="near", xaxt='n', xlab="", ylab="", col='blue')
    apply(cbind(nc, rollapply(kn$centers, width=5, mean, partial=T)), 1, function(x) { lines(x[-1], col=x[1], lwd=2) })
    axis(1, at=(1:ncol(near_n))[c(T,F,F)], labels=sapply(near_dists, n2str, 1)[c(T,F,F)], cex=0.7, las=2)

    plot(colMeans(far_n), type='l', lwd=4, ylim=range(kf$centers), main="far", xaxt='n', xlab="", ylab="", col='blue')
    apply(cbind(fc, kf$centers), 1, function(x) { lines(x[-1], col=x[1], lwd=2) })
    axis(1, at=(1:ncol(far_n))[c(T,F)], labels=sapply(far_dists, n2str, 1)[c(T,F)], cex=0.7, las=2)

    dev.off()
  }

  # grouping to phases
  mv$group1 = 3
  mv$group1[ mv$mitotic_bin > sch_mitotic_cutoff] = 1

  if (criteria == "orig") {
    mv$group1[ mv$group1 == 3 & mv$dom_score <  dom_cutoff & mv$far_tightness >= far_tight_cutoff] = 2
    mv$group1[ mv$group1 == 3 & mv$dom_score >= dom_cutoff & mv$high_near_f >= near_f_cuttoff] = 4
  }
  else if (criteria == "new") {
    mv$group1[ mv$group1 == 3 & mv$high_near_f <= sch_high_near_s_start] = 2
    mv$group1[ mv$group1 == 3 & mv$high_near_f >= sch_high_near_mid_s] = 4
  }
  else if (criteria == "newer") {
    mv$group1 = 4
    mv$group1[mv$f_near_band > sch_near_s_start & mv$f_near_band <= sch_near_mid_s] = 3
    mv$group1[mv$f_near_band <= sch_near_s_start] = 2
    mv$group1[mv$f_mitotic_band >= sch_mitotic_min_post_m & mv$f_near_band <= sch_near_max_post_m] = 1
    mv$group1[mv$f_near_band > sch_near_max_post_m  &  sch_pre_m_slope * mv$f_mitotic_band - mv$f_near_band + sch_pre_m_intercept < 0] = 5
    
  }


  # fix groups by clustering and reassigning cells based on cluster identity
  dec_nv = dec_n[rownames(mv), ]
  tfn = sprintf("/tmp/stam%s", Sys.getpid())

  ka = TGLKMeans_wrapper(dec_nv, tfn, round(nrow(dec_nv)/10))
  
  mv$kac = ka$cluster
  gr_t = t(table(mv$kac, mv$group1))

  #####gr_t is missing a cluster
  mv$group = rownames(gr_t)[apply(gr_t, 2, which.max)[as.character(mv$kac)]]
  #m[m$group1 == 1, 'group'] = 1

  if (plot_basic) {
    print(table(mv$group, mv$group1))
    if (criteria == "newer") {
      png(sprintf("%s/ordering/new_decay_groups.png", sch_fig_dir), 280*4, 280)
      par(mfrow=c(1,4))
      plot(mv$f_mitotic_band, mv$f_near_band, pch=19, cex=0.7, col=addalpha(mini_qualitative_colors[mv$group1], 0.6), xlab="%mitotic", ylab="%near", main="by criteria")
      grid(col='black', lwd=0.5)
      plot(mv$f_mitotic_band, mv$f_near_band, pch=19, cex=0.7, col=addalpha(mini_qualitative_colors[as.numeric(mv$group)], 0.6), xlab="%mitotic", ylab="%near", main="k-means fixed")
      grid(col='black', lwd=0.5)
      plot(mv$early_f/min(mv$early_f), mv$f_near_band,  pch=19, cex=0.7, col=addalpha(mini_qualitative_colors[as.numeric(mv$group)], 0.6), ylab="%near", xlab="repli-score", main="")
      grid(col='black', lwd=0.5)

      plot(colMeans(dec_nv), type='l', lwd=4, ylim=c(0,0.03), xlab="", ylab="", xaxt='n')
      li = seq(1, ncol(dec_nv), by=20)
      axis(1, at=li, labels=sapply(dec_dists, n2str, 1)[li], cex=0.7, las=1)
      for (gr in 1:max(mv$group)) {
        lines(colMeans(dec_nv[mv$group == gr,]), lwd=2, col=mini_qualitative_colors[gr])
      }
      grid(col='black')
      
      legend("topleft", legend=paste(c("post-M (", "G1 (", "early-S (", "S/G2 (", "pre-M"), table(mv$group), ")", sep=""), lwd=2, col=mini_qualitative_colors[1:5], bty='n', cex=0.8)

      dev.off()
    }
    else if (criteria == "new") {

      png(sprintf("%s/ordering/decay_by_phases.png", sch_fig_dir), 720, 280)
      layout(matrix(1:3, 1, 3))
      plot(colMeans(dec_nv), type='l', lwd=4, ylim=c(0,0.03), xlab="", ylab="", xaxt='n')
      li = seq(1, ncol(dec_nv), by=20)
      axis(1, at=li, labels=sapply(dec_dists, n2str, 1)[li], cex=0.7, las=1)
      lines(colMeans(dec_nv[mv$group == 1,]), lwd=2, col='orange')
      lines(colMeans(dec_nv[mv$group == 2,]), lwd=2, col='blue')
      lines(colMeans(dec_nv[mv$group == 3,]), lwd=2, col='green')
      lines(colMeans(dec_nv[mv$group == 4,]), lwd=2, col='red')
      grid(col='black')
      
      legend("topleft", legend=paste(c("M (", "G1 (", "transition (", "S ("), table(mv$group), ")", sep=""), lwd=2, col=c('orange', 'blue', 'green', 'red'), bty='n', cex=0.8)
      
      plot(mv$dom_score, mv$far_tightness, pch=19, cex=0.6, col=c("green", "red")[1 + (mv$near_f > near_f_cuttoff )], xlab="dom_score", ylab="far_tightness", main="orig")

      legend("topright", legend=paste("near ", c(">", "<"), near_f_cuttoff), pch=19, col=c("red", "green"), bty='n')
      abline(v=dom_cutoff, lty=2, lwd=2)
      segments(x0=0, x1=dom_cutoff, y0=far_tight_cutoff, y1=far_tight_cutoff, lwd=2, lty=2)
      
      plot(mv$high_near_f, mv$early_f, pch=19, cex=0.6, xlab="high_near_f", ylab="early_f", main="new")

      abline(v=c(sch_high_near_s_start, sch_high_near_mid_s), lty=2, lwd=2)

      dev.off()
    }
  }

  # scell decay heatmaps by groups + sorting
  #g1_ord = order(m[m$group == 1, 'far_wm'] + 1e-8 * m[m$group == 1, 'far_tightness'])

  #m_ord = order(m[m$group == 1, 'high_near_mu'])
  gr_ords = list()
  if (criteria == "orig") {
    gr_ords[[1]] = order(mv[mv$group == 1, 'all_mu'])
    gr_ords[[2]] = order(mv[mv$group == 2, 'far_mu'])
    gr_ords[[3]] = order(mv[mv$group == 3, 'high_near_f'] )
    gr_ords[[4]]  = order(mv[mv$group == 4, 'dom_score'])
  }
  else if (criteria == "new") {
    gr_ords[[1]] = order(mv[mv$group == 1, 'all_mu'])
    gr_ords[[2]] = order(mv[mv$group == 2, 'far_mu'])
    m3 = mv[mv$group == 3,]
    m4 = mv[mv$group == 4,]

    gr_ords[[3]] = order( (m3$high_near_f / var(m3$high_near_f)) + (m3$early_f / var(m3$early_f)))
    gr_ords[[4]]  = order( (m4$high_near_f / var(m4$high_near_f)) - (m4$early_f / var(m4$early_f)))
  }
  else if (criteria == "newer") {
    gr_ords[[1]] = order(mv[mv$group == 1, 'f_mitotic_band'], decreasing=T)
    gr_ords[[2]] = order(scale(mv[mv$group == 2, 'f_near_band']) + scale(mv[mv$group == 2, 'far_mu']))
    #repli_s = m$early_f / min(m$early_f)
    #repli_s = max(repli_s) - repli_s
    #gr_ords[[3]] = order(scale(m[m$group == 3, 'f_near_band']) - scale(repli_s[m$group == 3]))
    #gr_ords[[4]]  = order(scale(m[m$group == 4, 'f_near_band']) + scale(repli_s[m$group == 4]))
    near3 = mv[mv$group == 3, 'f_near_band']
    near4 = mv[mv$group == 4, 'f_near_band']
    tor3 = mv[mv$group == 3, 'early_f']
    tor4 = mv[mv$group == 4, 'early_f']
    gr_ords[[3]] = order(near3/var(near3) + tor3/var(tor3))
    gr_ords[[4]] = order(near4/var(near4) - tor4/var(tor4))

    gr_ords[[5]] = order(mv[mv$group == 5, 'f_mitotic_band'], decreasing=F)
    
  }

  all_cells = c()
  for (i in seq_along(gr_ords)) {
    all_cells = c(all_cells, rownames(mv)[mv$group  == i][ gr_ords[[i]] ])
  }

  cord = data.frame(x=seq_along(all_cells))
  rownames(cord) = all_cells


  mv$ord = cord[rownames(mv), 'x']

  all_cells = all_cells[mv[all_cells, 'valid']]
  
  message(nrow(mv), " valid cells of ", nrow(m))

  if (plot_extra) {
    op = plot_decay_mat_sorted_by_x(zlim=c(0, 0.025), order_by="given", name=sprintf("%s_valid_clus_", criteria), cells_names_ord=all_cells, remove_trans=T, add_dist_grid=T, add_ab=add_ab, add_nodig_trans=add_nodig_trans, remove_unear_dist=sch_remove_near_dist, vlines=cumsum(table(mv[all_cells, 'group']))/length(all_cells))


    ## png(sprintf("%s/ordering/dec_corr.png", sch_fig_dir), 2*nrow(mv)+50, 2*nrow(mv)+50)
    ## par(mar=c(1,1,1,1))
    ## dd = as.numeric(gsub("X", "", colnames(dec)))
    ## dd = dd[-length(dd)]
    ## dec_n = dec[, -c(1:last_to_remove, min(which(dd >= 1e8)):ncol(dec))]
    ## dec_n = dec_n / rowSums(dec_n)
    ## dnc = cor(t(dec_n[rownames(mv),]))
    ## diag(dnc) = NA

    ## #cor_col_spec = colorRampPalette(c(rep("blue", 6), "black", "yellow", "yellow", "red", "purple", "brown", "black"))(600)
    ## cor_col_spec = colorRampPalette(rev(c('black', '#a50026','#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695', '#313695')))(1200)
    ## image(dnc, zlim=c(-1, 1), col=cor_col_spec, main="", xaxt='n', yaxt='n', xlab="", ylab="")
    ## abline(v=cumsum(table(mv$group))/nrow(mv), lwd=3, col='green')
    ## abline(h=cumsum(table(mv$group))/nrow(mv), lwd=3, col='green')
    ## dev.off()

    ## plot_legend(sprintf("%s/ordering/dec_corr_leg.png", sch_fig_dir), c(-1,1), cor_col_spec, "pearson")
  }
      
  mv[all_cells, ]
}




#######
# insulation on TAD borders by CC order, grouped and border stratification (AA, BB, AB)
insu_by_cc_grouped <- function(cc_ord_cells_m, rel_odir_pref="ordering/insu_cc", scale=3e5, group_size=100, delta=1e-3, min_mean_b_count=1, chroms=paste0("chr", sch_chroms), idata=NULL, sort_domains="kmeans", n_ds=8, k=4, tor_horiz=300e3, norm_domains=T, norm_by_slices=NULL, colspec=c("cyan", "blue", "black", "yellow", "red"), name="", insu_name=sprintf("_tad_borders_gt%s", n2str(sch_remove_near_dist)), per_cell_func=mean, rebuild=F, dixon_max_border_size=80e3)
{
  cc_ord_cells_m = cc_ord_cells_m[cc_ord_cells_m$valid,]
  cc_ord_cells = rownames(cc_ord_cells_m)

  if (is.null(idata)) {
    bi = get_ds_insu_stats(scale=scale, res=0, coords_nm=insu_name, chroms=chroms, n_ds=n_ds, rebuild=rebuild)

    a = bi[['a']][, cc_ord_cells]
    b1 = bi[['b1']][, cc_ord_cells]
    b2 = bi[['b2']][, cc_ord_cells]

    ag = rollapply(t(a), width=group_size, by=group_size, FUN=sum)
    b1g = rollapply(t(b1), width=group_size, by=group_size, FUN=sum)
    b2g = rollapply(t(b2), width=group_size, by=group_size, FUN=sum)

    valid = colMeans(b1g) >= min_mean_b_count * group_size & colMeans(b2g) >= min_mean_b_count * group_size
    #valid = colMeans(ag + b1g + b2g) >= min_mean_b_count

    x_r = (ag[,valid] + delta) / (ag[,valid] + b1g[,valid] + b2g[,valid] + delta)
  }
  else {
    x_r = idata$x_r
    a = idata$a
    b1 = idata$b1
    b2 = idata$b2
  }

  if (norm_domains) {
    if (is.null(norm_by_slices)) {
      norm_by_slices = c(1, nrow(x_r))
    }
    #x  = t(log2( t(x_r)/ colMeans(x_r[norm_by_slices[1]:norm_by_slices[2],])))
    x = log2(x_r)
    x = t(t(x) - colMeans(x))
    #x = log2(t(t(x_r) / colMeans(x_r)))
    #colspec = rev(colspec)
  }
  else {
    x = log2(x_r)
  }

  sort_domains_str = switch(sort_domains,
    hclust="hclust",
    bord_tor_type="bord_tor_type",
    tor="tor",
    kmeans=paste0("kmeans_k", k))

  odir = sprintf("%s/%s%s_scale%s_grp%d_min%s_dom_by_%s%s", sch_fig_dir, rel_odir_pref, name, n2str(scale, 0), group_size, as.character(min_mean_b_count), sort_domains_str, ifelse(norm_domains, sprintf("_normDom%d-%d", norm_by_slices[1], norm_by_slices[2]), ""))
  system(paste("mkdir -p", odir))

  qs = quantile(x, c(0.5, 0.999))
  zlim = c( qs[1] - abs(diff(qs)), qs[1] + abs(diff(qs)))
  print(zlim)

  if (FALSE) { # old: classify TAD end instead of mid of boundary
    all_d = NULL
    d = .sch_get_pool_domains()
    d$early = ifelse(d$tor > 0, 1, 0)
    ds = split(d, as.character(d$chrom))
    for (chrom in names(ds)) {
      dc = ds[[chrom]]
      dc$up_early   = c(NA, d[1:(nrow(dc)-1), 'early'])
      dc$up_dist = c(1e8, d[2:nrow(dc), 'start'] - d[1:(nrow(dc)-1), 'end'] )
      dc$down_early = c(d[2:nrow(dc), 'early'], NA)
      dc$down_dist = c(d[2:nrow(dc), 'start'] - d[1:(nrow(dc)-1), 'end'], 1e8 )
      all_d = rbind(all_d, dc)
    }

    all_d$type_s = NA
    up_ind = 1:nrow(all_d) %% 2 == 1 & all_d$up_dist <= inter_dom_max_dist
    all_d[up_ind, 'type_s'] = all_d[up_ind, 'up_early'] + all_d[up_ind, 'early']

    all_d$type_e = NA
    down_ind = 1:nrow(all_d) %% 2 == 0 & all_d$down_dist <= inter_dom_max_dist
    all_d[down_ind, 'type_e'] = all_d[down_ind, 'down_early'] + all_d[down_ind, 'early']

    bord_tor = data.frame(tor=rep(all_d$tor, 2), chrom=rep(as.character(all_d$chrom), 2), coord=c(all_d$start, all_d$end), type=c(all_d$type_s, all_d$type_e))
    rownames(bord_tor) = c(paste0(as.character(d$chrom), "_X", d$start), paste0(as.character(d$chrom), "_X", d$end))
  }

  gvtrack.create("up_tor", sch_tor_track)
  gvtrack.iterator("up_tor", sshift=-tor_horiz, eshift=0)
  gvtrack.create("down_tor", sch_tor_track)
  gvtrack.iterator("down_tor", sshift=0, eshift=tor_horiz)

  if (length(grep("dixon", insu_name, ignore.case=T)) > 0) {
      b = read.delim(dixon_borders_fn, header=T)
      b = b[b$end - b$start <= dixon_max_border_size, ]
      b = intervals.centers(b)
    }
  else {
    b = .sch_get_pool_tad_borders()
  }

  btor = gextract("up_tor", "down_tor", gintervals.all(), iterator=b)
  bord_tor = data.frame(chrom=as.character(b$chrom), coord=b$start, type=ifelse(btor$up_tor > 0, 1, 0) + ifelse(btor$down_tor > 0, 1, 0), tor=(btor$up+btor$down)/2, up_tor=btor$up, down_tor=btor$down)
  rownames(bord_tor) = gsub("+", ".", paste0(as.character(b$chrom), "_X", b$start), fixed=T)

  for (chrom in chroms) {
    message("processing ", chrom)
    xc = x[, grep(paste0(chrom, "_X"), colnames(x))]
    xct = NULL

    png(sprintf("%s/ins_by_cc_%s.png", odir, chrom), height=ncol(xc)+500, 500)
    layout(matrix(2:1, 2, 1), height=c(1,4))
    par(mar=c(1,4,1,1))
    par(xaxs='i')

    if (sort_domains == "hclust") {
      dhc = hclust(dist(t(xc)), method="ward.D2")
      dom_ord = dhc$order
    }
    else if (sort_domains == "tor") {
      dom_ord = order(bord_tor[colnames(xc), 'tor'])
    }
    else if (sort_domains == "kmeans") {
      km = TGLKMeans_wrapper(t(xc), "/tmp/stam", k)
      dom_ord = order(km$cluster)
    }
    else if (sort_domains == "bord_tor_type") {
      dom_ord = order(bord_tor[colnames(xc), 'type'], na.last=F)
    }

    image(pmax(pmin(as.matrix(xc[, dom_ord]), zlim[2]), zlim[1]), zlim=zlim, col=colorRampPalette(colspec)(200), xaxt='n', yaxt='n', xlab="", ylab="")

    if (sort_domains == "kmeans") {
      blines = cumsum(km$size)/ncol(xc)
      abline(h=blines, col='red', lwd=4)
      axis(2, at=(c(0, blines[-k]) + blines)/2, labels=1:k, cex.axis=1)
      xct = km$cluster
    }
    else if (sort_domains == "bord_tor_type") {
      v = bord_tor[colnames(xc)[dom_ord], 'type']
      abline(h=cumsum(table(v))/ncol(xc), col='red', lwd=4)
      xct = bord_tor[colnames(xc), 'type']
    }

    plot(rowMeans(xc), type='l', lwd=2, xaxt='n', xlab="", ylab="", ylim=c(-4.5, -2))

    if (!is.null(xct)) {
      for (i in 1:max(xct)) {
        if (sum(xct == i) > 1) {
          v = rowMeans(xc[,xct == i])
        }
        else {
          v = xc[, xct == i]
        }
        lines(v, col=mini_qualitative_colors[i], lwd=2)
      }
      legend("bottomright", legend=1:max(xct), col=mini_qualitative_colors[1:max(xct)], lwd=2, lty=1, ncol=2, bty='n', cex=0.8)
    }

    dev.off()

    if (sort_domains == "kmeans") {
          png(sprintf("%s/ins_by_cc_%s_clus_tor.png", odir, chrom), height=400, 400)
          tab = table(km$cluster, bord_tor[names(km$cluster), 'type'])
          tab_n = t( tab / rowSums(tab))
          btype_nm = list("0"="LL", "1"="EL", "2"="EE")
          rownames(tab_n) = unlist(btype_nm[rownames(tab_n)])
          barplot(tab_n, legend.text=T, ylim=c(0,1.4))
          dev.off()
    }
  }

  # all data
  trends = NULL
  km = list()
  if (sort_domains == "kmeans") {
    km = TGLKMeans_wrapper(t(x), "/tmp/stam", k)
    xct = km$cluster
    dom_ord = order(km$cluster)
    }
  else if (sort_domains == "tor") {
    dom_ord = order(bord_tor[colnames(x), 'tor'])
    xct = NULL
  }
  else if (sort_domains == "bord_tor_type") {
    xct = bord_tor[colnames(x), 'type']
    km$cluster = xct
    names(km$cluster) = colnames(x)
    km$size = table(xct)

    dom_ord = order(xct, na.last=F)
  }
  ks = sort(unique(km$cluster))


  if (!is.null(xct)) {
    png(sprintf("%s/ins_by_cc_all_trends.png", odir), width=500, 400)
    trends = NULL
    #message(sprintf("clustered borders: %d of %d", length(xct), ncol(x)))
    for (i in ks) {
      if (sum(xct == i) > 1) {
        v = rowMeans(x[,xct == i])
      }
      else {
        v = x[, xct == i]
      }
      trends = rbind(trends, v)
    }
    rownames(trends) = paste0("c", ks)

    plot(rowMeans(x), type='l', lwd=4, xaxt='n', xlab="CC phase", ylab=ifelse(norm_domains, "insulation (norm)", "insulation"), ylim=range(trends), col=NA)

    apply(cbind(mini_qualitative_colors[1:nrow(trends)], trends), 1, function(y) { lines(y[-1], col=y[1], lwd=2) })

    kms = table(xct)
    legend("bottomright", legend=paste0(rownames(trends), "(", kms[ks], ")"), col=mini_qualitative_colors[1:nrow(trends)], lwd=2, lty=1, ncol=2, bty='n', cex=0.8)
    grid(nx=NA, ny=NULL, col='black')
    dev.off()
  }


  png(sprintf("%s/ins_by_cc_all_hm.png", odir), height=ncol(x)+100, width=800)
  image(pmax(pmin(as.matrix(x[, dom_ord]), zlim[2]), zlim[1]), zlim=zlim, col=colorRampPalette(colspec)(200), xaxt='n', yaxt='n', xlab="", ylab="")

  if (!is.null(xct)) {
    blines = cumsum(km$size)/ncol(x)
    abline(h=blines, col='red', lwd=4)
    axis(2, at=(c(0, blines[-nrow(trends)]) + blines)/2, labels=rownames(trends), cex.axis=3)
  }
  dev.off()

  if (!is.null(xct)) {
    png(sprintf("%s/ins_by_cc_all_clus_tor.png", odir), height=400, width=k*60)

    bord_tor$cluster = 0
    bord_tor[names(km$cluster), 'cluster'] = km$cluster

    ind = bord_tor$cluster >= 0
    tab = table(bord_tor$cluster[ind], bord_tor$type[ind])

    tab_n = t( tab / rowSums(tab))
    btype_nm = list("0"="LL", "1"="EL", "2"="EE")
    rownames(tab_n) = unlist(btype_nm[rownames(tab_n)])
    barplot(tab_n, legend.text=T, ylim=c(0,1.4))
    dev.off()
  }

  # insulation per cluster x cell by CC order

  #layout(matrix(seq_along(ks), length(ks), 1))
  yc = NULL
  mean_ins_per_cell = log2(apply(a + delta, 2, per_cell_func)/apply(a + b1 + b2 + delta, 2, per_cell_func))
  ind = split(names(km$cluster), km$cluster)

  for (i in seq_along(ks)) {
    yc = rbind(yc, log2(apply(a[ind[[i]], ] + delta, 2, per_cell_func)/apply(a[ind[[i]], ] + b1[ind[[i]], ] + b2[ind[[i]], ] + delta, 2, per_cell_func)))
  }

  for (i in seq_along(ks)) {

    png(sprintf("%s/ins_per_clus_and_cell_by_cc_c%s.png", odir, ks[i]), height=400, width=ncol(a)/3+200)
    par(mar=c(1,3,1,1))

    plot(mean_ins_per_cell, pch=19, cex=0.5, main="", xlab="", ylab="A / all (log2)", ylim=range(yc))
    lines(cyclic_rollmean(mean_ins_per_cell, width=51), lwd=4)

    points(yc[i,], pch=19, cex=0.5, col=mini_qualitative_colors[i])
    lines(cyclic_rollmean(yc[i,], width=51), lwd=4, col=mini_qualitative_colors[i])
    legend("topleft", legend=paste0("c", ks[i]), lwd=4, lty=1, col=mini_qualitative_colors[i])
    dev.off()
  }

  cw = 101
  ycm = cyclic_rollmean(yc, width=cw)
  png(sprintf("%s/ins_per_clus_and_cell_by_cc_trends.png", odir), height=400, width=ncol(a)/3+200)
  #layout(matrix(1:2, 2, 1))
  par(mar=c(1,3,1,1))

  plot(ycm[1,], type='l', lty=1, lwd=4, main="", xlab="", ylab="A / all (log2)", ylim=range(ycm), col=mini_qualitative_colors[1])
  apply(cbind(mini_qualitative_colors[2:nrow(ycm)], ycm[-1, ]), 1, function(v) { lines(v[-1], lwd=4, col=v[1]) })
  grid(col='black')
  legend("bottomright", legend=paste0("c", ks), lwd=4, lty=1, col=mini_qualitative_colors[seq_along(ks)], cex=2)
  dev.off()

  ycr = -t(log2(t(yc) / colMeans(yc)))

  ycr_m = cyclic_rollmean(ycr, width=cw)
  ycr_sd = cyclic_rollapply(ycr, width=cw, func_str="sd")
  ycr_n = cyclic_rollapply(ycr, width=cw, func_str="length")
  ycr_se = ycr_sd / sqrt(ycr_n)
  n_sd = 1.5
  outliers = which(colSums(ycr >= ycr_m + n_sd * ycr_sd | ycr <= ycr_m - n_sd * ycr_sd) >= 0.75*nrow(ycr))

  message(length(outliers), " outliers")
  png(sprintf("%s/ins_per_clus_and_cell_by_cc_mean_rel_trends.png", odir), height=nrow(ycr)*220, width=ncol(a)/3+200)
  layout(matrix(1:nrow(ycr), nrow(ycr), 1))
  par(mar=c(1,3,1,1))

  for (i in 1:nrow(ycr)) {
    plot(ycr_m[i,], type='l', lwd=4, main="", xlab="", ylab="A / all (log2)", ylim=range(ycr[i,]), col=mini_qualitative_colors[i])
    points(ycr[i,], pch=19, cex=0.2)
    lines(ycr_m[i,] + n_sd * ycr_sd[i,], lwd=2, lty=2, col=mini_qualitative_colors[i])
    lines(ycr_m[i,] - n_sd * ycr_sd[i,], lwd=2, lty=2, col=mini_qualitative_colors[i])
    grid(col='black')
    #abline(v=outliers, col='red', lty=2)
    legend("bottomright", legend=paste0("c", ks[i]), lty=1, lwd=4, col=mini_qualitative_colors[i], cex=1)
  }
  dev.off()

  png(sprintf("%s/ins_per_clus_and_cell_by_cc.png", odir), height=400, width=ncol(a)/3+200)
  par(mar=c(1,3,1,1))
  plot(mean_ins_per_cell, pch=19, cex=0.7, main="", xlab="", ylab="A / all (log2)", ylim=range(yc))
  rect(xleft=(1:length(mean_ins_per_cell))-0.5, xright=(1:length(mean_ins_per_cell))+0.5, ybottom=apply(yc, 2, min), ytop=apply(yc, 2, max), border=NA, col='darkgray')
  lines(cyclic_rollmean(mean_ins_per_cell, width=51), lwd=4)
  points(mean_ins_per_cell, cex=0.7, pch=19)
  grid(col='gray')
  dev.off()

  png(sprintf("%s/ins_per_clus_and_cell_by_cc_cvar.png", odir), height=400, width=ncol(a)/3+200)
  par(mar=c(1,3,1,1))
  y = apply(yc, 2, max) - apply(yc, 2, min)
  plot(y, pch=19, cex=0.7, main="", xlab="", ylab="intra-clus diff")
  lines(cyclic_rollmean(y, width=51), lwd=4)
  grid(col='gray')
  dev.off()

  png(sprintf("%s/ins_per_clus_and_cell_by_cc_distr.png", odir), height=400, width=400)
  plot(density(y), main="", xlab="cell border insulation", ylab="density", lwd=0, ylim=range(apply(yc, 1, function(v) { range(density(v)$y) } )), xlim=range(yc))
  sapply(1:nrow(yc), function(s) { lines(density(yc[s,]), col=mini_qualitative_colors[as.numeric(s)], lwd=2) })
  legend("topleft", legend=paste0("c", ks), lty=1, lwd=2, col=mini_qualitative_colors[1:nrow(yc)])
  dev.off()

  if (sort_domains == "kmeans") {
    cg = sort(unique(cc_ord_cells_m$group))
    gc = split(rownames(cc_ord_cells_m[colnames(yc), ]), cc_ord_cells_m[colnames(yc), 'group'])
    combs = combn(x=nrow(yc), m=2)
    for (i in 1:ncol(combs)) {
      c1 = combs[1, i]
      c2 = combs[2, i]

      nr = floor(sqrt(length(cg)))
      nc = ceiling(length(cg) / nr)

      png(sprintf("%s/ins_per_clus_and_cell_c%s_c%s.png", odir, ks[c1], ks[c2]), height=220*nr, width=220*nc)
      layout(matrix(1:nrow(yc), nr, nc, byrow=T))
      par(mar=c(4,4,1,1))

      for (g in cg) {
        plot(yc[c1,], yc[c2,], pch=19, col='gray', cex=0.5, xlim=range(yc), ylim=range(yc), main="", xlab=paste0("c", ks[c1]), ylab=paste0("c", ks[c2]))

        points(yc[c1, gc[[g]]], yc[c2, gc[[g]]], pch=19, col=colorRampPalette(c("orange", "red"))(length(gc[[g]])), cex=0.5)
        grid(col='black')
        abline(a=0, b=1, lty=3)
        legend("topleft", legend=paste0("g", g), pch=19, cex=0.5, col="red")
      }

      dev.off()

    }
  }

  list(x_r=x_r, x=x, bord_tor=bord_tor, trends=trends, dom_ord=dom_ord, a=a, b1=b1, b2=b2, km=km, yc=yc, outliers=colnames(ycr)[outliers], mean_ins_per_cell=mean_ins_per_cell)
}





classify_contacts_by_ctcf <- function(contacts, ctcf_tns=sch_ctcf_tns, extend_by=5e3, ctcf_rep_q=0.999, chroms=sch_chroms, motif_threshold=13, tor_horiz=200e3, ctcf=NULL)
{

  .classify_side = function(v) {
    v = unique(v)
    ifelse ("B" %in% v || "F" %in% v && "R" %in% v, "B", ifelse("F" %in% v, "F", ifelse("R" %in% v, "R", ifelse("C" %in% v, "C", "N"))))
  }

  if (is.null(ctcf)) {
    message("finding ctcf intervals")
    ctcf_vtns = paste0("v_", ctcf_tns)
    for (i in seq_along(ctcf_tns)) {
      gvtrack.create(ctcf_vtns[i], ctcf_tns[i], "global.percentile.max")
    }

    chip = gscreen(paste(paste(ctcf_vtns, ">", ctcf_rep_q), collapse=" | "), intervals=gintervals(chroms), iterator=50)

    gvtrack.create("ctcf_e_forward", "motifs.CTCF_forward", "global.percentile.max")
    gvtrack.create("ctcf_e_reverse", "motifs.CTCF_reverse", "global.percentile.max")
                                        #ctcf = gextract("log10(ctcf_e_forward)", "log10(ctcf_e_reverse)", chip, iterator=50, colnames=c("forward", "reverse"))

    gvtrack.create("up_tor", sch_tor_track)
    gvtrack.iterator("up_tor", sshift=-tor_horiz, eshift=0)
    gvtrack.create("down_tor", sch_tor_track)
    gvtrack.iterator("down_tor", sshift=0, eshift=tor_horiz)

    ctcf = gextract("-log2(1 - ctcf_e_forward)", "-log2(1 - ctcf_e_reverse)", "up_tor", "down_tor", chip, iterator=50, colnames=c("forward", "reverse", "up_tor", "down_tor"))

    ctcf$pos_m = ctcf$forward >= motif_threshold
    ctcf$neg_m = ctcf$reverse >= motif_threshold

    ctcf$type = "C"
    ctcf[ctcf$pos_m, 'type'] = 'F'
    ctcf[ctcf$neg_m, 'type'] = 'R'
    ctcf[ctcf$pos_m & ctcf$neg_m, 'type'] = "B"
  }

  if (!is.null(contacts)) {
    c1 = data.frame(chrom=contacts$chrom1, start=contacts$start1, end=contacts$end1, id=1:nrow(contacts))
    c2 = data.frame(chrom=contacts$chrom2, start=contacts$start2, end=contacts$end2, id=1:nrow(contacts))

    n1 = gintervals.neighbors(c1, ctcf, maxneighbors=1e6, mindist=-extend_by, maxdist=extend_by, na.if.notfound=T)
    n2 = gintervals.neighbors(c2, ctcf, maxneighbors=1e6, mindist=-extend_by, maxdist=extend_by, na.if.notfound=T)

    contacts$type1 = tapply(n1$type, n1$id, .classify_side)
    contacts$type2 = tapply(n2$type, n2$id, .classify_side)

    contacts$type = "one_sided"

    contacts[contacts$type1 == 'R' & contacts$type2 == 'R' | contacts$type1 == 'F' & contacts$type2 == 'F', 'type'] = "same"
    contacts[contacts$type1 == 'R' & contacts$type2 == 'F' |
             contacts$type1 == 'B' & contacts$type2 == 'F' |
             contacts$type1 == 'R' & contacts$type2 == 'B' |
             contacts$type1 == 'B' & contacts$type2 == 'B', 'type'] = "conv"
    contacts[contacts$type1 == 'F' & contacts$type2 == 'R' |
             contacts$type1 == 'B' & contacts$type2 == 'R' |
             contacts$type1 == 'F' & contacts$type2 == 'B', 'type'] = "diver"
    contacts[contacts$type1 == 'C' | contacts$type2 == 'C', 'type'] = "chip"
    contacts[contacts$type1 == 'N' & contacts$type2 == 'N', 'type'] = "none"
  }

  list(ctcf=ctcf, contacts=contacts)

}

#########
# contacts between promoters (stratified by gene activity)
.sch_classify_tss <- function(rna_tns=c("rna.129_Cast.ES.pf_rep1", "rna.129_Cast.ES.pf_rep2"), bps_into_gene=1000, enr_cutoff=6, k4me3_around_tss=2e3, extend_tss_for_union=5e3)
{
  for (i in seq_along(rna_tns)) {
    gvtrack.create(paste0("rna_gpm", i), rna_tns[i], "global.percentile.max")
    gvtrack.iterator(paste0("rna_gpm", i), sshift=0, eshift=bps_into_gene)
  }

  tss = gintervals.load("intervs.global.tss")
  tss = tss[ grep("chrM|chrY", as.character(tss$chrom), perl=T, invert=T), ]

  tss$rna = gextract(sprintf("-log2(1-pmax(%s))", paste0("rna_gpm", seq_along(rna_tns), collapse=",")), intervals=tss, iterator=tss, colnames="x")$x

  tss$type = ifelse(tss$rna >= enr_cutoff, "T", "F")

  gvtrack.create("k4me3_r1",  "Encode.esb4.h3k4me3.rep1", "global.percentile.max")
  gvtrack.create("k4me3_r2",  "Encode.esb4.h3k4me3.rep2", "global.percentile.max")

  gvtrack.iterator("k4me3_r1", sshift=-k4me3_around_tss, eshift=k4me3_around_tss)
  gvtrack.iterator("k4me3_r2", sshift=-k4me3_around_tss, eshift=k4me3_around_tss)

  tss$k4me3 = gextract("-log2(1 - pmax(k4me3_r1, k4me3_r2))", intervals=tss, iterator=tss, colnames="x")$x

  ctss = gintervals.canonic(gintervals.force_range(data.frame(chrom=tss$chrom, start=tss$start - extend_tss_for_union, end=tss$end + extend_tss_for_union)))
  tss$map = attr(ctss, "mapping")

  tss$id = 1:nrow(tss)

  ind = do.call(c, lapply(split(tss, tss$map), function(x) { x[which.max(x$k4me3)[1], 'id']}))

  tss$selected = F
  tss[ind, 'selected'] = T

  tss
}


############
# features should be 1d interval, with a 'type' column that will be used for stratification
plot_contacts_around_feature <- function(tns, features, features_2d=NULL, ofn="ttt.png", width_per_panel=180, extend_by=25e3, chroms = sch_chroms, band = c(2e5, 1e6), shift_x_features_by=0, shift_y_features_by=0, name="", glob_zlim=NULL, raw_z=NULL, norm_by_z=NULL, colspec=c("blue", "white", "red"), smooth_method="bins", nbins=51, raw_counts=T, x_types=NULL, y_types=NULL, rebuild=F, mar=c(3,3,2,2), plot_axis=T, plot_inline=F)
{
  stopifnot(is.null(features) | is.null(features_2d))
  
  if (!is.null(features_2d)) {
    features = data.frame(
      chrom= c(as.character(features_2d$chrom1), as.character(features_2d$chrom2)),
      start= c(features_2d$start1, features_2d$start2),
      end=   c(features_2d$end1, features_2d$end2),
      strand=c(features_2d$strand1, features_2d$strand2),
      type=  c(features_2d$type1, features_2d$type2))

  }

  types = sort(unique(features$type))
  if (is.null(x_types)) {
    x_types = types
  }
  if (is.null(y_types)) {
    y_types = types
  }

  raw_z_fn = sprintf("%s/c_around_f_%d_tns_%d_f_exby_%s_chrs_%s_band_%s-%s-shift%d.%d_name_%s_%s_%d_bins_x%s_y%s.Rdata", sch_rdata_dir, length(tns), nrow(features), n2str(extend_by, 0), paste0(chroms, collapse="-"), n2str(band[1], 0), n2str(band[2], 0), shift_x_features_by, shift_y_features_by, name, smooth_method, nbins, paste0(sort(x_types), collapse=""), paste0(sort(y_types), collapse=""))


  if (is.null(raw_z)) {
    if (!file.exists(raw_z_fn) | rebuild) {
      features_x = features[features$type %in% x_types, ]
      if (shift_x_features_by > 0) {
        message("shifting X features upstream by ", shift_x_features_by)
        features_x$start = pmax(features_x$start - shift_x_features_by, 0)
        features_x$end   = features_x$end - shift_x_features_by
        features_x = features_x[features_x$start < features_x$end, ]
      }

      features_y = features[features$type %in% y_types, ]
      if (shift_y_features_by > 0) {
        message("shifting Y features upstream by ", shift_y_features_by)
        features_y$start = pmax(features_y$start - shift_y_features_by, 0)
        features_y$end   = features_y$end - shift_y_features_by
        features_y = features_y[features_y$start < features_y$end, ]
      }


      mat_filler = expand.grid(1:nbins, 1:nbins)
      colnames(mat_filler) = c('bin1', 'bin2')

      dd = NULL
      for (i in seq_along(tns)) {
        tn = tns[i]
        message(i, "\t", tn)
        x = gextract(tn, gintervals.2d(chroms), band=-rev(band), colnames='cov')
        x1 = data.frame(chrom=x$chrom1, start=x$start1, end=x$end1, nn_id=1:nrow(x))
        x2 = data.frame(chrom=x$chrom2, start=x$start2, end=x$end2, nn_id=1:nrow(x))
        t1 = gintervals.neighbors(x1, features_x, maxneighbors=1, mindist=-extend_by, maxdist=extend_by)
        t2 = gintervals.neighbors(x2, features_y, maxneighbors=1, mindist=-extend_by, maxdist=extend_by)
        ids = as.character(intersect(t1$nn_id, t2$nn_id))
        rownames(t1) = as.character(t1$nn_id)
        rownames(t2) = as.character(t2$nn_id)
        if (!is.null(features_2d)) {
          t12 = cbind(t1[ids, 5:6], t2[ids, 5:6], ids)
          colnames(t12) = c('chrom1', 'start1', 'chrom2', 'start2', 'ids')
          t12_filt = merge(t12, features_2d, by=c('chrom1', 'start1', 'chrom2', 'start2'))
          ids = as.character(t12_filt$ids)
        }

        dd = rbind(dd, cbind(t1[ids, c('dist', 'type')], t2[ids, c('dist', 'type')]))
      }
      colnames(dd) = c('dist1', 'type1', 'dist2', 'type2')


      binsize = 2 * extend_by / nbins
      dd$bin1 = pmin(ceiling(pmax(dd$dist1 + extend_by, 1) / binsize), nbins)
      dd$bin2 = pmin(ceiling(pmax(dd$dist2 + extend_by, 1) / binsize), nbins)
      raw_z = list()

      for (x in x_types) {                        
        for (y in y_types) {
          ind = dd$type1 == x & dd$type2 == y
          message(sprintf("%s%s: %d contacts (out of %d)", x, y, sum(ind), nrow(dd)))
          if (smooth_method == "kde") {
            kdd = kde2d(dd[ind, 'dist1'], dd[ind, 'dist2'], n=nbins)
            raw_z[['x']] = kdd$x
            raw_z[['y']] = kdd$y
            z = kdd$z
          }
          else if (smooth_method == "bins") {
            bins = seq(-extend_by, extend_by, length=nbins)
            raw_z[[paste0(x,y)]] = table(rbind(mat_filler, dd[ind, c('bin1', 'bin2')])) - 1
            raw_z[['x']] = bins
            raw_z[['y']] = bins
          }
        }
      }
      save(raw_z, file=raw_z_fn)
    }
  }
  else {
    save(raw_z, file=raw_z_fn)
  }
  load(raw_z_fn)
  #message(sprintf("loaded %s", raw_z_fn))

  nx = length(x_types)
  ny = length(y_types)

  if (!is.null(ofn) | plot_inline) {
    if (!plot_inline) {
      fig_png(ofn, nx * width_per_panel, ny*width_per_panel)
      par(mfrow=c(nx, ny))
      par(mar=mar)
    }


    for (x in x_types) {
      for (y in y_types) {
        z = raw_z[[paste0(x, y)]]
        if (sum(z) > 0) {
          zn = z
          if (raw_counts) {
            zn = z / sum(z)
          }

          if (!is.null(norm_by_z)) {
            norm_z = norm_by_z[[paste0(x, y)]]
            #message(sprintf("norm mat has %d contacts", sum(norm_z)))
            if (raw_counts) {
              norm_z = norm_z + 1
              zn = z + 1
              zn = log2((zn / sum(zn)) / (norm_z / sum(norm_z)))
            }
            else {
              zn = log2(zn / norm_z)
            }
          }

          #message(sprintf("data range for %s%s: %g to %g\tq(0.01, 0.99): %g to %g", x, y, min(zn), max(zn), quantile(zn, 0.01), quantile(zn, 0.99)))

          if (is.null(glob_zlim)) {
            zlim = range(zn)
          }
          else {
            zlim = glob_zlim
          }


          image(raw_z$x, raw_z$y, pmin(pmax(zn, zlim[1]), zlim[2]), zlim=zlim, col=colorRampPalette(colspec)(1000), main=ifelse(plot_axis, sprintf("%s %s%s %s(%g-%g)", name, x, y, ifelse(raw_counts, sprintf("(%d) ", sum(z)), ""), zlim[1], zlim[2]), ""),  xlab="", ylab="", xaxt=ifelse(plot_axis, 's', 'n'), yaxt=ifelse(plot_axis, 's', 'n'))
          #grid(col='black')
        }
        else {
          plot(1:10, col=NA, main=sprintf("%s %s%s", name, x, y), xaxt='n', yaxt='n', xlab="", ylab="")
          text(x=1, y=5, labels="NO DATA", pos=4)
        }


      }
    }

    if (!plot_inline) {
      dev.off()
    }
  }

  raw_z
}


#######
# mitotic alignment (in trans contacts)
sch_mit_align_test <- function(rebuild=F, nms=NULL)
{
  if (is.null(nms)) { 
    mv = sch_decay_metrics
    nms = rownames(mv)
  }

  fn = sprintf("%s/trans_contacts_corr.tab", sch_table_dir)

  if (!file.exists(fn) || rebuild) {
    align = NULL

    commands = paste0(".sch_align_nm(nm=\"", nms, "\")", collapse=",")

    res = eval(parse(text=paste("gcluster.run(",commands,")")))

    for (i in seq_along(res)) {
      rv = res[[i]]$retv
      stopifnot(!is.null(rv), typeof(rv) == "double")
      align = rbind(align, rv)
    }

    rownames(align) = nms
    align = cbind(align, nms)
    colnames(align) = c('n', 'ntrans', 'c_align', 'frac_pos', 'nm')

    write.table(align, fn, quote=F, sep="\t")
  }
  align = read.table(fn, header=T, stringsAsFactors=F)

  sw=101
  png(sprintf("%s/trans_corr_by_cell_sw%d.png", sch_fig_dir, sw), 960, 350)
  plot(align$c_align, pch=19, cex=0.5, ylab="corr", main="", xlab="")
  lines(cyclic_rollmean(as.numeric(align$c_align), width=sw), lwd=4)
  #abline(v=cumsum(table(m$group)), col='blue', lwd=2, lty=2)
  grid(col='black', nx=NA, ny=NULL)
  dev.off()
}

.sch_align_nm <- function(nm)
{
  chrs = gintervals.all()
  rownames(chrs) = chrs$chrom

  a = gextract(nm, gintervals.2d.all())

  ainter = a[as.character(a$chrom1) != as.character(a$chrom2),]

  ainter$start1 = ainter$start1/chrs[as.character(ainter$chrom1), "end"]
  ainter$start2 = ainter$start2/chrs[as.character(ainter$chrom2), "end"]

  xs = split(ainter, paste(ainter$chrom1, ainter$chrom2))
  xc = unlist(lapply(xs, function(m) { if (nrow(m) > 10) { cor(m$start1, m$start2) } else { NA }  } ))

  c(nrow(a), nrow(ainter), cor(ainter$start1, ainter$start2), mean(xc > 0, na.rm=T))
}

#####################
plot_contacts_between_score_filtered_conv_ctcf_by_pooled_group <- function(extend_by=50e3, nbins=51, ctcf_chip_q=0.9, tor_horiz=2e5, band=c(2e5, 1e6), pool_q_cutoff=0, pool_loop_extend=1e4, min_score=60, score_tn=sprintf("%s_score", pool_tn), conv_ctcf=NULL, zlim=c(-1,1), colspec=c("purple", "navy", "blue", "#87FFFF", "white", "#FF413D", "black", "orange", "yellow"), odir=sprintf("%s/ordering", sch_fig_dir))
{
  group_nums = 1:5
  group_nms = c('post_m', 'g1', 'early_s', 'mid_s_g2', 'pre_m')
  
  group_tns = c(paste(pool_tn, "group", group_nums, group_nms, sep="_"))

  if (is.null(conv_ctcf)) {
    conv_ctcf = .get_conv_ctcf_with_tor(ctcf_chip_q=ctcf_chip_q, tor_horiz=tor_horiz, band=band, pool_q_cutoff=pool_q_cutoff, pool_loop_extend=pool_loop_extend, min_d_score=min_score, score_tn=score_tn)

  }
  conv_ctcf$strand1 = conv_ctcf$strand2 = 1
  conv_ctcf$type1 = conv_ctcf$type2 = "A"


  ps_ctcf_RvsF_dd = plot_contacts_around_feature(sprintf("%s_shuffle", pool_tn), features=NULL, features_2d=conv_ctcf, extend_by=extend_by, nbins=nbins, ofn=sprintf("%s/pool_shuffled_conv_ctcf_%dbins_min%d_obs_%s.png", odir, nbins, min_score, n2str(extend_by, 0)), name=paste0("ps_conv_ctcf_m", min_score), width_per_panel=300)

  for (i in seq_along(group_nums)) {
    g_dd = plot_contacts_around_feature(group_tns[i], features=NULL,
      features_2d=conv_ctcf, extend_by=extend_by, nbins=nbins,
      ofn=sprintf("%s/g%s_cells_conv_ctcf_%dbins_min%d_norm_%s.png", odir, group_nums[i], nbins, min_score, n2str(extend_by, 0)), width_per_panel=300,
      name=sprintf("min%d_g%d_%s", min_score, group_nums[i], group_nms[i]),
      norm_by_z=ps_ctcf_RvsF_dd, colspec=colspec, glob_zlim=zlim)
  }

  plot_legend(sprintf("%s/gAll_cells_conv_ctcf_%dbins_min%d_norm_%s.png", odir, nbins, min_score, n2str(extend_by, 0)), zlim=zlim, crp=colorRampPalette(colspec)(200), name="enr")

  conv_ctcf
}


######################
#
gen_domains_decay_by_size_and_tor <- function(dom_dist_cutoffs=c(0, 17.9, 18.9, 19.9, 20.5), decay_step=0.1, chroms=sch_chroms, dom_cutoff_q=0.15, rebuild=F, per_len_bin_ylim=T, ylim_qs=c(0.005, 0.995), width_per_panel=450, height_per_panel=130, tor_cols=list(EE="red", E="pink", L="lightblue", LL="blue"))
{

  tor_types = c("EE", "E", "L", "LL")
  rdata_ofn = sprintf("%s/intra_domain_decay_%s_by_%.2f_chroms_%s_%f.Rdata", sch_rdata_dir, paste(dom_dist_cutoffs, collapse="-"), decay_step, paste(chroms, collapse="-"), dom_cutoff_q)

  if (file.exists(rdata_ofn) && !rebuild) {
    load(rdata_ofn)
  }
  else {
    cells = rownames(sch_decay_metrics)[sch_decay_metrics$valid]

    d = .sch_get_pool_domains(dom_cutoff_q=dom_cutoff_q)
    d = d[!is.na(d$tor), ]

    d$tor_bin = "EE"
    d$tor_bin[d$tor > 0 & d$tor <= median(d$tor[d$tor > 0])] = "E"
    d$tor_bin[d$tor > median(d$tor[d$tor <= 0]) & d$tor <= 0] = "L"
    d$tor_bin[d$tor <= median(d$tor[d$tor <= 0])] = "LL"

    d = d[d$len < max(2**dom_dist_cutoffs), ]
    for (i in 2:length(dom_dist_cutoffs)) {
      d[d$len >= 2**dom_dist_cutoffs[i-1] & d$len < 2**dom_dist_cutoffs[i], 'len_bin'] = i-1
    }

    slen = sapply(2**dom_dist_cutoffs, n2str, 2)
    d$len_sbin = paste(slen[-length(slen)], slen[-1], sep="_")[d$len_bin]
    write.table(table(d$tor_bin, d$len_sbin), sprintf("%s/intra_domain_decay_l%s_q%f_domain_counts.tab", sch_table_dir, paste(dom_dist_cutoffs, collapse="-"), dom_cutoff_q), sep="\t", quote=F)

    intra = list()
    for (len_s in 1:max(d$len_bin)) {
      intra[[len_s]] = list()

      decay_breaks = 2**seq(log2(sch_remove_near_dist), dom_dist_cutoffs[len_s + 1], by=decay_step)
      message(sprintf("decay breaks: %f (%f) to %f (%f)", min(decay_breaks), log2(min(decay_breaks)), max(decay_breaks), log2(max(decay_breaks))))

      for (tor_s in tor_types) {

        c_intra = NULL
        c_d = d[d$len_bin == len_s & d$tor_bin == tor_s, ]

        commands = paste0("gcis_decay(\"", cells, "\", decay_breaks, src=gintervals(chroms), domain=c_d, intervals=gintervals.2d(chroms), include.lowest=F)", collapse=",")

        res = eval(parse(text=paste("gcluster.run(",commands,")")))

        for (i in seq_along(res)) {
          rv = res[[i]]$retv
          stopifnot(!is.null(rv), typeof(rv) == "double")
          c_intra = rbind(c_intra, rv[,1])
        }

        rownames(c_intra) = cells
        colnames(c_intra) = log2(decay_breaks[-length(decay_breaks)])
        print(range(c_intra))

        intra[[len_s]][[tor_s]] = c_intra
      }
    }
    save(intra, file=rdata_ofn)
  }

  for (i in 1:length(intra)) {
    df = NULL
    for (tor_s in tor_types) {
      x = intra[[i]][[tor_s]]
      x_n = x / rowSums(x)
      x_n[rowSums(x) == 0, ] = 0
      v = as.matrix(x_n) %*% 2^as.numeric(colnames(x_n))
      df = rbind(df, as.vector(v))
    }
    rownames(df) = tor_types
    if (per_len_bin_ylim) {
      glob_ylim = quantile(df, ylim_qs)
    }
    else {
      glob_ylim = NULL
    }
    .plot_points_and_trend(t(df), sprintf("%s/ordering/intra_domain_decay_by_cc_d%d_%s-%s%s.png", sch_fig_dir, i, n2str(2**dom_dist_cutoffs[i], 1), n2str(2**dom_dist_cutoffs[i+1], 1), ifelse(per_len_bin_ylim, "_gYlim", "")), cols=NULL,  height_per_panel=height_per_panel, width_per_panel=width_per_panel, sw=101, col_by_cond=F, collapse_q=ylim_qs, cex=0.4, glob_ylim=glob_ylim, lcol="blue", pcol="black", plot_trends=T, trend_cols=unlist(tor_cols[tor_types]))
  }

  intra
}


###############
#
chrom_phasing_stats_by_cc <- function(chroms=sch_chroms, pointsize=7, fres=120, sw=51, odir=sprintf("%s/ordering", sch_fig_dir), ex_chroms=c(1, 11), cex=0.5, width=240, height_per_panel=90, pcol='red', trend_col='black', cex_lab=1.4)

{

  mv = sch_decay_metrics[sch_decay_metrics$valid, ]

  x = as.matrix(sch_chrom_decay_res[rownames(mv), ])

  x = x[, grep("trans", colnames(x), invert=T)]
  x_dists = as.numeric(gsub("chr.*_", "", colnames(x), perl=T))
  x = x[, x_dists >= sch_remove_near_dist ]

  far_mu = NULL
  near_f = NULL

  for (chr in chroms) {
    cind = grep(paste0("chr", chr, "_"), colnames(x))
    xc = x[, cind]

    colnames(xc) = gsub("chr.*_", "", colnames(xc), perl=T)

    xc_dists = as.numeric(colnames(xc))

    far_cols = xc_dists >= sch_far_dists[1] & xc_dists <= sch_far_dists[2]
    far = xc[, far_cols]
    far_dists = xc_dists[far_cols]
    far_n = far / rowSums(far)

    far_mu = rbind(far_mu, t(as.matrix(far_n) %*% far_dists))

    near_cols = xc_dists >= sch_near_dists[1] & xc_dists <= sch_near_dists[2]
    near = xc[, near_cols]
    near_dists = xc_dists[near_cols]

    near_f = rbind(near_f, rowSums(near) / rowSums(xc))
  }

  rownames(far_mu) = rownames(near_f) = chroms
  colnames(far_mu) = colnames(near_f) = rownames(mv)

  near_f = near_f[, rownames(mv)[mv$group %in% 2:4]]
  far_mu = far_mu[, rownames(mv)[mv$group == 2]]

  chr_ints = gintervals(chroms)
  chr_lens = chr_ints$end
  names(chr_lens) = chr_ints$chrom

  near_f_o = near_f[order(chr_lens[paste0('chr', rownames(near_f))]), ]

  far_mu_l = log2(far_mu)
  far_mu_l_o = far_mu_l[order(chr_lens[paste0('chr', rownames(near_f))]), ]


  near_f_o_s   = t(rollapply(t(near_f_o), width=sw, partial=T, FUN=mean))
  far_mu_l_o_s = t(rollapply(t(far_mu_l_o), width=sw, partial=T, FUN=mean))

  png(sprintf("%s/chrom_near_frac_by_cc_sw%d.png", odir, sw), width=width, height=height_per_panel * (1 + length(ex_chroms)), pointsize=pointsize, res=fres)
  layout(1:(length(ex_chroms)+1))
  par(mar=c(2,2,0,1))
  ylim = range(near_f[as.character(ex_chroms), ])
  xseq = seq(min(mv[colnames(near_f), 'ord']), by=1, len=ncol(near_f))

  plot(xseq, colMeans(near_f), type='l', ylim=ylim, col=NA, xlab="", ylab="%near", xaxt='n', cex.lab=cex_lab)
  apply(cbind(colorRampPalette(c("lightgray", "black"))(nrow(near_f_o_s)), near_f_o_s), 1, function(v) { lines(xseq, v[-1], col=v[1], lwd=1) } )
  grid(col='black', lwd=0.5)

  for (chr in ex_chroms) {
    plot(xseq, near_f[as.character(chr), ], ylim=ylim, pch=19, cex=cex, col=pcol, xlab="", ylab="%near", main="", xaxt='n', cex.lab=cex_lab)
    lines(xseq, near_f_o_s[as.character(chr), ], lwd=2, col=trend_col)
    if (chr == ex_chroms[length(ex_chroms)]) {
      axis(1)
    }
    grid(col='black', lwd=0.5)
    legend("bottomright", legend=paste0("chr", chr), col=NA, bty='n', cex=1.6)
  }
  dev.off()

  png(sprintf("%s/chrom_far_mean_contact_dist_by_cc_sw%d.png", odir, sw), width=width, height=height_per_panel * (1 + length(ex_chroms)), pointsize=pointsize, res=fres)
  layout(1:(length(ex_chroms)+1))
  par(mar=c(2,2,0,1))
  ylim = range(far_mu_l[as.character(ex_chroms), ])
  xseq = seq(min(mv[colnames(far_mu), 'ord']), by=1, len=ncol(far_mu))
  plot(xseq, colMeans(far_mu_l_o), type='l', ylim=ylim, col=NA, xlab="", ylab="avg far dist", xaxt='n', cex.lab=cex_lab)
  apply(cbind(colorRampPalette(c("lightgray", "black"))(nrow(far_mu_l_o_s)), far_mu_l_o_s), 1, function(v) { lines(xseq, v[-1], col=v[1], lwd=1) } )
  grid(col='black', lwd=0.5)

  for (chr in ex_chroms) {
    plot(xseq, far_mu_l[as.character(chr), ], ylim=ylim, pch=19, cex=cex, col=pcol, xlab="", ylab="avg far dist", main="", xaxt='n', cex.lab=cex_lab)
    lines(xseq, far_mu_l_o_s[as.character(chr), ], lwd=2, col=trend_col)
    if (chr == ex_chroms[length(ex_chroms)]) {
      axis(1)
    }
    grid(col='black', lwd=0.5)
    legend("bottomright", legend=paste0("chr", chr), col=NA, bty='n', cex=1.6)
  }
  dev.off()

  list(far_mu=far_mu, near_f=near_f)
}


##########################################
# conv CTCF contacts at a band, strarify by CTCF ToR
sch_conv_ctcf_contact_count_by_tor <- function(tor_horiz=2e5, ctcf_chip_q=0.99, loop_width=10e3, band=c(200e3, 1e6), conv_ctcf=NULL, nm_from=1, nm_size=NULL, job_mode=F, chroms=sch_chroms, conv_ctcf_pool_q=0, min_d_score=60, count_contact_once=T, name="")
{
  if (job_mode) {
    options(gmultitasking=F)
    options(gmax.data.size=1e9)
  }

  if (is.null(conv_ctcf)) {
    message("recalc conv ctcf loops...")
    conv_ctcf = .get_conv_ctcf_with_tor(ctcf_chip_q, tor_horiz, band, conv_ctcf_pool_q, loop_width, min_d_score=min_d_score)
  }

  res = NULL

  nms = rownames(sch_decay_metrics)[sch_decay_metrics$valid] #sch_decay_metrics)[sch_decay_metrics$valid]

  nm_to = ifelse(is.null(nm_size), length(nms), min(nm_from+nm_size-1, length(nms)))
  nms = nms[nm_from:nm_to]

  #commands = paste0(".sch_inner_contacts_by_ftype(\",", nms, "\", conv_ctcf, loop_width, 'tor')", collapse=",")

  #res = eval(parse(text=paste("gcluster.run(",commands,")")))

  m = matrix(0, length(nms), 16)
  rownames(m) = nms
  types =  c('EE', 'EL', 'LE', 'LL')
  colnames(m) = apply(as.matrix(expand.grid(types, types)), 1, paste0, collapse="")

  regions = c('total', 'oo', 'oi', 'ii', 'io')
  m_l = list()
  for (reg in regions) {
    m_l[[reg]] = m
  }

  for (i in seq_along(nms)) {
    message(sprintf("processing %d (%s) on %d loops", i, nms[i], nrow(conv_ctcf)))
    r = .sch_inner_contacts_by_ftype(nms[i], conv_ctcf, loop_width, 'tor', chroms, band, count_contact_once)
    for (reg in regions) {
      m_l[[reg]][nms[i], names(r[[reg]])] = r[[reg]]
    }
    #r = res[[i]]
    #if (typeof(r$retv) == "integer") {
    #  m[nms[i], names(r$retv)] = r$retv
    #}
  }

  for (reg in regions) {
    ofn = sprintf("%s/%sconv_ctcf_min%d_%s_counts_%d_to_%d_loop_size%s_q%s%s.tab", sch_rdata_dir, name, min_d_score, reg, nm_from, nm_to, n2str(loop_width, 0), n2str(conv_ctcf_pool_q), ifelse(count_contact_once, "_one2one", "_one2many"))
    message(paste("writing to", ofn))
    write.table(m_l[[reg]], ofn, quote=F, sep="\t")
  }

  list(ctcf=conv_ctcf, res=m_l)
}

####
.get_conv_ctcf_with_tor <- function(ctcf_chip_q=0.99, tor_horiz=2e5, band=c(2e5, 1e6), pool_q_cutoff=0, pool_loop_extend=1e4, pool_loop_extend_bg=3e4, pool_loop_foc_enr_q=0.8, min_d_score=60, score_tn=paste0(pool_tn, "_score"), only_non_overlap_loops=T, external_loops=NULL)
{
  if (is.null(external_loops)) {
    ct = classify_contacts_by_ctcf(NULL, ctcf_rep_q=ctcf_chip_q, tor_horiz=tor_horiz)
    ctcf = ct$ctcf[ct$ctcf$type != 'C' & !is.na(ct$ctcf$up_tor) & !is.na(ct$ctcf$down_tor) ,]

    ctcf$tor = paste0(ifelse(ctcf$up_tor > 0, 'E', 'L'), ifelse(ctcf$down_tor > 0, 'E', 'L'))

    conv_ctcf = gintervals.neighbors(ctcf[ctcf$type == 'R' | ctcf$type == 'B',], ctcf[ctcf$type == 'F' | ctcf$type == 'B',], maxneighbors=1e9, mindist=band[1], maxdist=band[2])
    colnames(conv_ctcf) = c(paste0(colnames(ctcf), '1'), paste0(colnames(ctcf), '2'), 'dist')
    conv_ctcf = conv_ctcf[, c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'tor1', 'tor2')]
    
    # remove divergent loops
    conv_ctcf = filter(conv_ctcf, start1 < start2)
  }
  
  else {
    # loops supplied, adding tor data to anchors
    conv_ctcf = intervals.2d.centers(external_loops)
    anchors = unique(data.frame(chrom=rep(as.character(conv_ctcf$chrom1), 2), start=c(conv_ctcf$start1, conv_ctcf$start2), end=c(conv_ctcf$end1, conv_ctcf$end2)))

    gvtrack.create("up_tor", sch_tor_track)
    gvtrack.iterator("up_tor", sshift=-tor_horiz, eshift=0)
    gvtrack.create("down_tor", sch_tor_track)
    gvtrack.iterator("down_tor", sshift=0, eshift=tor_horiz)
    
    tor = gextract("up_tor", "down_tor", intervals=anchors, iterator=anchors, colnames=c('up', 'down'))
    tor = tor[!is.na(tor$up) & !is.na(tor$down), ]
    
    tor$tor = paste0(ifelse(tor$up > 0, 'E', 'L'), ifelse(tor$down > 0, 'E', 'L'))
   

    
    conv_ctcf = merge(conv_ctcf, tor[, c('chrom', 'start', 'tor')], by.x=c('chrom1', 'start1'), by.y=c('chrom', 'start'))
    conv_ctcf = merge(conv_ctcf, tor[, c('chrom', 'start', 'tor')], by.x=c('chrom1', 'start2'), by.y=c('chrom', 'start'), suffixes=c('1', '2'))
    conv_ctcf = conv_ctcf[, c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'tor1', 'tor2')]

    loop_dist = conv_ctcf$start2 - conv_ctcf$start1
    conv_ctcf = conv_ctcf[loop_dist >= band[1] & loop_dist <= band[2], ]
  }
    
  
  if (!is.na(min_d_score) && is.null(external_loops)) {
    message(sprintf("filtering %d loops by pool track score (>= %d)", nrow(conv_ctcf), min_d_score))
    gvtrack.create("s_tn", score_tn, "max")
    gvtrack.iterator.2d("s_tn", -pool_loop_extend, pool_loop_extend, -pool_loop_extend, pool_loop_extend)
    conv_ctcf$score = gextract("ifelse(is.na(s_tn), -100, s_tn)", iterator=conv_ctcf, intervals=conv_ctcf, colnames='score')$score
    #print(quantile(conv_ctcf$score, seq(0, 1, by=0.005)))
    conv_ctcf = conv_ctcf[conv_ctcf$score >= min_d_score, ]

  }

  if (pool_loop_extend > 0 && pool_loop_extend_bg > 0 && is.null(external_loops)) {
    message(sprintf("filtering %d loops by pool track loop foci enrichment (%f -ile)", nrow(conv_ctcf), pool_loop_foc_enr_q))
    gvtrack.create("foc", pool_tn, "weighted.sum")
    gvtrack.create("bg", pool_tn, "weighted.sum")
    gvtrack.iterator.2d("foc", -pool_loop_extend, pool_loop_extend, -pool_loop_extend, pool_loop_extend)
    gvtrack.iterator.2d("bg", -pool_loop_extend_bg, pool_loop_extend_bg, -pool_loop_extend_bg, pool_loop_extend_bg)
    foc = gextract("ifelse(is.na(foc), 0, foc)", "ifelse(is.na(bg), 0, bg)", iterator=conv_ctcf, intervals=conv_ctcf, colnames=c('foc', 'bg'))
    conv_ctcf$enr = log2( (foc$foc + 1) / (foc$bg + 1))
    conv_ctcf = conv_ctcf[conv_ctcf$enr >= quantile(conv_ctcf$enr, pool_loop_foc_enr_q), ]
  }
  
  if (pool_q_cutoff > 0 && is.null(external_loops)) {
    message(sprintf("filtering %d loops by pool track coverage (%f -ile)", nrow(conv_ctcf), pool_q_cutoff))
    gvtrack.create("ptws", pool_tn, "weighted.sum")
    gvtrack.iterator.2d("ptws", -pool_loop_extend, pool_loop_extend, -pool_loop_extend, pool_loop_extend)

    cov = gextract("ifelse(is.na(ptws), 0, ptws)", intervals=gintervals.2d.all(), iterator=conv_ctcf, colnames="count")

    conv_ctcf$cov = cov$count
    #conv_ctcf$enr = log2(cov$count / cov$back_count)
    q_min = quantile(conv_ctcf$cov, pool_q_cutoff)
    #q_min_enr = quantile(conv_ctcf$enr, pool_q_cutoff)
    message(sprintf("Filtering %d loops to min %d contacts", nrow(conv_ctcf), q_min))
    conv_ctcf = conv_ctcf[conv_ctcf$cov >= q_min, ]

  }

  if (only_non_overlap_loops && is.null(external_loops)) {
    message(sprintf("filtering %d loops to non-overlapping ones (at %s)", nrow(conv_ctcf), n2str(pool_loop_extend_bg, 0)))

    nn = gintervals.neighbors.wrap(conv_ctcf, conv_ctcf, maxneighbors=2)
    nn = nn[nn$start1_2 != nn$start1_1 | nn$start2_1 != nn$start2_2, ]
    nn_filt = nn[nn$dist1 >= pool_loop_extend_bg | nn$dist2 >= pool_loop_extend_bg, ]
    conv_ctcf = nn_filt[, 1:10]
    colnames(conv_ctcf) = gsub("_1", "", colnames(conv_ctcf))

  }
  message(sprintf("returning %d loops", nrow(conv_ctcf)))
          
  conv_ctcf
}

#####
.sch_inner_contacts_by_ftype <- function(x, f, maxdist=1e4, type_name='tor', chroms=sch_chroms, band=c(200e3, 1e6), count_contact_once=T)
{
  xc = gextract(x, intervals=gintervals.2d(chroms), band=c(-band[2]-maxdist, -band[1] + maxdist))

  nn = gintervals.neighbors.wrap(xc, f, maxneighbors=ifelse(count_contact_once, 1, 1e6), mindist1=-maxdist, maxdist1=maxdist, mindist2=-maxdist, maxdist2=maxdist)

  # o = out of TAD, i = in
  nn_oo = nn[nn$start1_1 - nn$start1_2 >= -maxdist & nn$start1_1 - nn$start1_2 <= 0 & nn$start2_1 - nn$start2_2 >= 0 & nn$start2_1 - nn$start2_2 <= maxdist , ]
  nn_oi = nn[nn$start1_1 - nn$start1_2 >= -maxdist & nn$start1_1 - nn$start1_2 <= 0 & nn$start2_1 - nn$start2_2 >= -maxdist & nn$start2_1 - nn$start2_2 <= 0 , ]
  nn_ii = nn[nn$start1_1 - nn$start1_2 >= 0 & nn$start1_1 - nn$start1_2 <= maxdist & nn$start2_1 - nn$start2_2 >= -maxdist & nn$start2_1 - nn$start2_2 <= 0 , ]
  nn_io = nn[nn$start1_1 - nn$start1_2 >= 0 & nn$start1_1 - nn$start1_2 <= maxdist & nn$start2_1 - nn$start2_2 >= 0 & nn$start2_1 - nn$start2_2 <= maxdist , ]

  cn1 = paste0(type_name, '1_2')
  cn2 = paste0(type_name, '2_2')

  total = oo = oi = ii = io = NULL
  if (!is.null(nn))    {
    total = table(paste0(nn[, cn1], nn[, cn2]))
  }
  if (!is.null(nn_oo)) {
    oo    = table(paste0(nn_oo[, cn1], nn_oo[, cn2]))
  }
  if (!is.null(nn_oi)) {
    oi    = table(paste0(nn_oi[, cn1], nn_oi[, cn2]))
  }
  if (!is.null(nn_ii)) {
    ii    = table(paste0(nn_ii[, cn1], nn_ii[, cn2]))
  }
  if (!is.null(nn_io)) {
    io    = table(paste0(nn_io[, cn1], nn_io[, cn2]))
  }

  res = list(total=total, oo=oo, oi=oi, ii=ii, io=io)
  lapply(res, function(x) { message(sum(x)) })
  res
}

##################
#
plot_conv_ctcf_counts <- function(width=410, height_per_panel=140, sw=101, cex=0.4, collapse_q=c(0.002, 0.998), col_by_cond=T, ctcf_chip_q=0.99, tor_horiz=2e5, band=c(2e5, 1e6), conv_ctcf=NULL, n_contact_scale=1e6, loop_widths=c(1e4, 3e4), conv_ctcf_pool_q=0, group_size_for_boxplot=150, boxplot_comp_inds=NULL, show_pvals_as_nums=T, fn_pref=sprintf("%s/ordering/", sch_fig_dir), rebuild=F, plot_horiz=F, min_d_score=60, score_tn=paste0(pool_tn, "_score_k100"), cells_per_job=10, en_boxplot_ylim=NULL, haploid=F, loop_foc_enr_q=0.8, count_contact_once=T, only_non_overlap_loops=T, external_loops=NULL, name="", lwd=2)
{

  if (is.null(conv_ctcf)) {

    conv_ctcf = .get_conv_ctcf_with_tor(ctcf_chip_q, tor_horiz, band, conv_ctcf_pool_q, pool_loop_extend=loop_widths[1], pool_loop_extend_bg=loop_widths[2], pool_loop_foc_enr_q=loop_foc_enr_q, min_d_score=min_d_score, score_tn=score_tn, only_non_overlap_loops=only_non_overlap_loops, external_loops=external_loops)

  }

  n_pairs = table(paste0(conv_ctcf$tor1, conv_ctcf$tor2))
  print(n_pairs)
  
  mv = sch_decay_metrics[sch_decay_metrics$valid, ]

  if (rebuild) {
    for (lw in loop_widths) {
      file.remove(list.files(path=sch_rdata_dir, pattern=sprintf("%sconv_ctcf_min%d_total_counts_.*_to_.*_loop_size%s_q%s.tab", name, min_d_score, n2str(lw, 0), n2str(conv_ctcf_pool_q)), full.names=T))

      nm_from = seq(1, nrow(mv), by=cells_per_job)
      commands = paste0("sch_conv_ctcf_contact_count_by_tor(tor_horiz=tor_horiz, ctcf_chip_q=ctcf_chip_q, loop_width=lw, band=band, conv_ctcf=conv_ctcf, nm_size=cells_per_job, job_mode=T, chroms=sch_chroms, conv_ctcf_pool_q=conv_ctcf_pool_q, min_d_score=min_d_score, nm_from=", nm_from, ", count_contact_once=count_contact_once, name=name)", collapse=",")
      res = eval(parse(text=paste("gcluster.run(",commands,")")))
    }
  }


  y = do.call("rbind", lapply(list.files(path=sch_rdata_dir, pattern=sprintf("%sconv_ctcf_min%d_total_counts_.*_to_.*_loop_size%s_q%s%s.tab", name, min_d_score, n2str(loop_widths[1], 0), n2str(conv_ctcf_pool_q), ifelse(count_contact_once, "_one2one", "_one2many")), full.names=T), read.table, header=T))


  yb = do.call("rbind", lapply(list.files(path=sch_rdata_dir, pattern=sprintf("%sconv_ctcf_min%d_total_counts_.*_to_.*_loop_size%s_q%s%s.tab", name, min_d_score, n2str(loop_widths[2], 0), n2str(conv_ctcf_pool_q), ifelse(count_contact_once, "_one2one", "_one2many")), full.names=T), read.table, header=T))

  message(sprintf("read from %sconv_ctcf_min%d_total_counts_.*_to_.*_loop_size%s_q%s%s.tab: mean EEEE: %f , %f",  name, min_d_score, n2str(loop_widths[1], 0), n2str(conv_ctcf_pool_q), ifelse(count_contact_once, "_one2one", "_one2many"), mean(y[,'EEEE']), mean(yb[,'EEEE'])))

  y = y[rownames(mv), ]
  yb = yb[rownames(mv), ]
  
  print(mean(as.matrix(y), na.rm=T))
  
  dn = sch_glob_decay_res[rownames(mv), ]
  dn = dn[, c(-1, -ncol(dn))]
  dn_dists = as.numeric(gsub("X", "", colnames(dn)))
  vind = dn_dists >= sch_remove_near_dist
  dn = dn[, vind]
  dn_dists = dn_dists[vind]

  dn_n = dn / rowSums(dn)

  cc = rowSums(dn)

  y_n = (n_contact_scale * y) / cc
  y_n = y_n[, names(n_pairs)]
  y_nn = t(t(y_n) / as.vector(n_pairs))

  x = data.frame("EE_EE"=y_nn[,"EEEE"], "LL_LL"=y_nn[,"LLLL"], "EE_LL"=rowSums(y_nn[, c("EELL", "LLEE")]), "EE_EL"=rowSums(y_nn[, c("EEEL", "LEEE")]), "EL_LL"=rowSums(y_nn[, c("ELLL", "LLLE")]))

  y = y + 1
  yb = yb + 1
  
  x_en = data.frame(
    "EE_EE"=y[,"EEEE"]/yb[,"EEEE"],
    "LL_LL"=y[,"LLLL"]/yb[,"LLLL"],
    "EE_LL"=rowSums(y[, c("EELL", "LLEE")])/rowSums(yb[, c("EELL", "LLEE")]),
    "EE_EL"=rowSums(y[, c("EEEL", "LEEE")])/rowSums(yb[, c("EEEL", "LEEE")]),
    "EL_LL"=rowSums(y[, c("ELLL", "LLLE")])/rowSums(yb[, c("ELLL", "LLLE")]))
  #x_en = log2(x_en)

  cols = get_cond_cols(haploid)

  .plot_points_and_trend(x, ofn=sprintf("%scontacts_in_conv_ctcf_%s-%s_min%d_by_tor_sw%d_q%s.png", fn_pref, n2str(band[1], 0), n2str(band[2], 0), min_d_score, sw, n2str(conv_ctcf_pool_q)), cols, width_per_panel=width, height_per_panel=height_per_panel, sw=sw, cex=cex, lwd=lwd, col_by_cond=col_by_cond, collapse_q=collapse_q, lcol='black')

  exp_enr = (loop_widths[1]/loop_widths[2])**2
  x_en_n = log2(x_en / exp_enr)
  x_en_n [ x_en_n < -1] = -1
  print(quantile(x_en_n, collapse_q, na.rm=T))

  for (show_pvals in c(T,F)) {
    .plot_points_and_trend(x_en_n, ofn=sprintf("%scontacts_in_conv_ctcf_%s-%s_min%d_by_tor_sw%d_%s_vs_%s_q%s%s_en.png", fn_pref, n2str(band[1], 0), n2str(band[2], 0), min_d_score, sw, n2str(loop_widths[1], 0), n2str(loop_widths[2], 0), n2str(conv_ctcf_pool_q), ifelse(count_contact_once, "_one2one", "_one2many")), cols, width_per_panel=width, height_per_panel=height_per_panel, sw=sw, cex=cex, lwd=lwd, col_by_cond=col_by_cond, collapse_q=collapse_q, lcol='black', glob_ylim=en_boxplot_ylim, group_size_for_boxplot=group_size_for_boxplot, selected_cols=c('EE_EE', 'LL_LL'), boxplot_comp_inds=boxplot_comp_inds, show_pvals_as_nums=show_pvals_as_nums, show_pvals=show_pvals, vertical=!plot_horiz, boxplot_grid=T)

   # .plot_points_and_trend(x_en, ofn=sprintf("%s/ordering/contacts_in_conv_ctcf_by_tor_sw%d_%s_vs_%s_q%s_en_constY.png", sch_fig_dir, sw, n2str(loop_widths[1], 0), n2str(loop_widths[2], 0), n2str(conv_ctcf_pool_q)), cols, width=width, height_per_panel=height_per_panel, sw=sw, col_by_cond=col_by_cond, collapse_q=collapse_q, lcol='black', glob_trend=rep(exp_enr, nrow(x_en)), glob_col='red', glob_ylim=c(exp_enr, quantile(unlist(x_en), collapse_q[2])), group_size_for_boxplot=group_size_for_boxplot, selected_cols=c('EE_EE', 'LL_LL'), boxplot_comp_inds=boxplot_comp_inds, show_pvals_as_nums=show_pvals_as_nums, show_pvals=show_pvals)
  }

  png(sprintf("%scontacts_in_conv_ctcf_%s-%s_min%d_by_tor_sw%d_q%s%s_trends.png", fn_pref, n2str(band[1], 0), n2str(band[2], 0), min_d_score, sw, n2str(conv_ctcf_pool_q), ifelse(count_contact_once, "_one2one", "_one2many")), 600, 300)

  x_s = cyclic_rollmean(t(x), sw)

  tcols = list("EE_EE"="red", "LL_LL"="blue", "EE_LL"="green", "EE_EL"="orange", "EL_LL"="cyan")

  plot(x_s[1, ], col=NA, ylim=range(x_s, na.rm=T), xlab='', ylab="norm contacts")
  apply(cbind(unlist(tcols[colnames(x)]), x_s), 1, function(v) { lines(v[-1], col=v[1], lwd=2) })
  grid(col='black', lwd=0.5)
  legend("topleft", legend=colnames(x), col=unlist(tcols[colnames(x)]), lty=1, lwd=2, bty='n', ncol=2)
  dev.off()

  # plot ii/oo enrichment
  oo = do.call("rbind", lapply(list.files(path=sch_rdata_dir, pattern=sprintf("conv_ctcf_min%d_oo_counts_.*_to_.*_loop_size%s_q%s.tab", min_d_score, n2str(loop_widths[2], 0), n2str(conv_ctcf_pool_q)), full.names=T), read.table, header=T))
  ii = do.call("rbind", lapply(list.files(path=sch_rdata_dir, pattern=sprintf("conv_ctcf_min%d_ii_counts_.*_to_.*_loop_size%s_q%s.tab", min_d_score, n2str(loop_widths[2], 0), n2str(conv_ctcf_pool_q)), full.names=T), read.table, header=T))
  oi = do.call("rbind", lapply(list.files(path=sch_rdata_dir, pattern=sprintf("conv_ctcf_min%d_oi_counts_.*_to_.*_loop_size%s_q%s.tab", min_d_score, n2str(loop_widths[2], 0), n2str(conv_ctcf_pool_q)), full.names=T), read.table, header=T))
  io = do.call("rbind", lapply(list.files(path=sch_rdata_dir, pattern=sprintf("conv_ctcf_min%d_io_counts_.*_to_.*_loop_size%s_q%s.tab", min_d_score, n2str(loop_widths[2], 0), n2str(conv_ctcf_pool_q)), full.names=T), read.table, header=T))
  x_quart = cbind((ii$EEEE+1e-3)/oo$EEEE, (ii$LLLL+1e-3)/oo$LLLL)
  colnames(x_quart) = c('EE', 'LL')

  #for (show_pvals in c(T, F)) {
  #  .plot_points_and_trend(x_quart, ofn=sprintf("%scontacts_in_conv_ctcf_min%d_by_tor_sw%d_q%s_%s_ii_vs_oo.png", fn_pref, min_d_score, sw, n2str(conv_ctcf_pool_q), n2str(loop_widths[2], 0)), cols, width_per_panel=width, height_per_panel=height_per_panel, sw=sw, col_by_cond=col_by_cond, collapse_q=collapse_q, lcol='blue', glob_trend=cyclic_rollmean(rowSums(ii)/rowSums(oo), sw), glob_col='black', group_size_for_boxplot=group_size_for_boxplot, boxplot_comp_inds=boxplot_comp_inds, show_pvals_as_nums=show_pvals_as_nums, show_pvals=show_pvals, vertical=!plot_horiz)
  #}

  list(x=x, x_enr=x_en_n, conv_ctcf=conv_ctcf, raw1=y, raw2=yb, n_pairs=n_pairs, oo=oo, ii=ii, oi=oi, io=io)
}



#===================================

sch_create_domain_marg_cov <- function(discard_cis_below=sch_remove_near_dist, count_dups=F, ds_by_exp_gcontent=F, n1_quant=0.1)
{
  mv = sch_decay_metrics[sch_decay_metrics$valid, ]
  mv$exp_gcontent = 1
  mv[mv$group == 1 | mv$group == 5 , 'exp_gcontent'] = 2
  s_ind = which(mv$group == 3 | mv$group == 4)
  mv[s_ind, 'exp_gcontent'] = seq(1, 2, length=length(s_ind))

  dn = sch_glob_decay_res[rownames(mv), ]
  dn_dists = as.numeric(gsub("X", "", colnames(dn)))
  mv$valid_cov = rowSums(dn[,  is.na(dn_dists) | dn_dists >= sch_remove_near_dist])

  n1_base_count = quantile(mv$valid_cov[mv$exp_gcontent == 1], n1_quant)
  mv$ds_to = round(n1_base_count * mv$exp_gcontent)

  if (ds_by_exp_gcontent) {
    nms = rownames(mv) [ mv$ds_to <= mv$valid_cov ]
    message(sprintf("downsamplign from %d cells out of %d", length(nms), nrow(mv)))
  }
  else {
    nms = rownames(mv)
    mv$ds_to = 0
  }
  mv = mv[nms, ]

  d = .sch_get_pool_domains()

  stats = matrix(0, length(nms), nrow(d))
  rownames(stats) = nms
  colnames(stats) = as.character(d$id)

  commands = paste0(".sch_gen_domain_marg_cov(nm=\"", nms, "\", ds_to=", mv$ds_to, sprintf(", discard_cis_below=%f, d=d, count_dups=%s)", discard_cis_below, ifelse(count_dups, "T", "F")), collapse=",")

  res = eval(parse(text=paste("gcluster.run(",commands,")")))

  for (i in 1:length(res)) {
    nm = nms[i]
    rv = res[[i]]$retv
    if (!is.null(rv)) {
      stopifnot(typeof(rv) == "integer")
      stats[nm, names(rv)] = rv
    }
  }

  if (ds_by_exp_gcontent) {
    sch_domain_raw_marg_cov_ds <<- stats
    write.table(sch_domain_raw_marg_cov_ds, sprintf("%s/domain_raw_marg_cov_ds%.2f.txt", sch_table_dir, n1_quant), quote=F, sep="\t")
  }
  else {
    sch_domain_raw_marg_cov <<- stats
    write.table(sch_domain_raw_marg_cov, sprintf("%s/domain_raw_marg_cov.txt", sch_table_dir), quote=F, sep="\t")
  }

}

####
.sch_gen_domain_marg_cov <- function(nm, ds_to=0,  discard_cis_below=sch_remove_near_dist, d=NULL, count_dups=F)
{
  if (is.null(d)) {
    d = .sch_get_pool_domains()
  }

  res = NULL

  x = gextract(nm, intervals=gintervals.2d.all(), colnames='count')
  x = x[as.character(x$chrom1) != as.character(x$chrom2) | abs(x$start1 - x$start2) >= sch_remove_near_dist, ]

  if (!is.null(x)) {
    if (ds_to > 0) {
      stopifnot(!count_dups)
      x = x[sample(1:nrow(x), ds_to), ]
    }

    x1 = data.frame(chrom=x$chrom1, start=x$start1, end=x$end1, count=x$count)

    xa = gintervals.neighbors.wrap(x1, d, maxneighbors=1, mindist=0, maxdist=0)

    if (!is.null(xa)) {
      res = tapply(xa$count_1, xa$id_2, ifelse(count_dups, sum, length))
    }
  }
  res
}


domain_marg_cov_plots <- function(mean_kb_cov_range=c(0.04, 0.16), clus_method="sort_by_cor", hc_method="ward.D2", order_by="tor", k=40, zlim=c(-0.3, 0.3), reord=NULL, chroms=sch_chroms, name="all", colspec=c("blue", "white", "red"), group_size=211, sw=101, max_chr_cov_enr=0.3, ord_only=F, odir=sch_fig_dir, plot_fig=T, use_method="raw", ord_pc_ifn=sprintf("%s/new_pc.txt", sch_table_dir), rebuild=F, ref_tn="scell.nextera.pool_good_hyb_2i_all_es")
{
  mv = sch_decay_metrics[sch_decay_metrics$valid, ]
  km = NULL
  d = .sch_get_pool_domains()

  use_cells = rownames(sch_chrom_marg_enr[ apply(sch_chrom_marg_enr, 1, function(v) max(abs(v))) <= max_chr_cov_enr, ])

  use_cells = intersect(use_cells, sch_good_cells)
  dm = sch_domain_raw_marg_cov[use_cells, ]
  #colnames(dm) = gsub("X", "", colnames(dm))

  stopifnot(sum(colnames(dm) == d$id) == ncol(dm))

  rownames(d) = d$id

  a_s = calc_domains_a_score(tn=ref_tn, rebuild=rebuild)

  if (order_by == "tor") {
    ord_vec =  -d$tor
  } else if (order_by == "trans_A_score") {
    ord_vec = -a_s[colnames(dm), 'trans_A_score']
  }
  else {
    ord_vec = nrow(d):1
  }
  names(ord_vec) = colnames(dm)

  ord_ind = !is.na(ord_vec)

  ord_vec = ord_vec[ord_ind]
  d = d[ord_ind, ]
  dm = dm[, ord_ind]
  a_s = a_s[ord_ind, ]
  
  ind = is.element(as.character(d$chrom), paste0('chr', chroms))
  ord_vec = ord_vec[ind]
  d = d[ind, ]
  dm = dm[, ind]
  a_s = a_s[ind, ]

  mpkb = colMeans(dm) / (d$len/1e3)
  message(sprintf("filtering %d of %d contacts (%.3f)", sum(mpkb >= mean_kb_cov_range[1] & mpkb <= mean_kb_cov_range[2]), ncol(dm), sum(mpkb >= mean_kb_cov_range[1] & mpkb <= mean_kb_cov_range[2])/ ncol(dm)))
  dm = dm[, mpkb >= mean_kb_cov_range[1] & mpkb <= mean_kb_cov_range[2]]
  
  dm = dm[, permute(1:ncol(dm))]

  d = d[colnames(dm), ]
  ord_vec = ord_vec[colnames(dm)]
  a_s = a_s[colnames(dm), ]

  dm_n = dm / rowSums(dm)

  dcor = cor(dm_n)

  if (file.exists(ord_pc_ifn)) {
    eord = read.table(ord_pc_ifn, header=T)
    eord$n = 1:nrow(eord)
    eord = eord[intersect(rownames(eord), colnames(dm_n)), ]
  }

  if (order_by == "tor" | order_by == "trans_A_score") {
    eg = dm_n[, tail(names(ord_vec)[order(ord_vec)], max(50, group_size))]
    lg = dm_n[, head(names(ord_vec)[order(ord_vec)], max(50, group_size))]
    message(sprintf("mean early %.2f\tmean late %.2f", mean(ord_vec[colnames(eg)], na.rm=T), mean(ord_vec[colnames(lg)], na.rm=T)))
  }
    

  if (ord_only) {
    reord = order(ord_vec)
  }

  if (is.null(reord)) {
    if (clus_method == "kmeans") {
      message("kmean clustering...")
      km = TGLKMeans_wrapper(dcor, sprintf("%s/dcor_kmeans", sch_rdata_dir), k)

      d$cl = km$cluster[as.character(d$id)]

      ec = which.max(tapply(d$tor, d$cl, mean, na.rm=T))
      lc = which.min(tapply(d$tor, d$cl, mean, na.rm=T))
      grps = lapply(1:k, function(cl) { dm_n[,  as.character(rownames(d)[!is.na(d$cl) & d$cl == cl])] })
      message(paste(ec, lc))
      v = tapply(d$tor, d$cl, mean, na.rm=T)
      ev = sapply(1:k, function(cl) { mean(cor(grps[[ec]], grps[[cl]])) })
      lv = sapply(1:k, function(cl) { mean(cor(grps[[lc]], grps[[cl]])) })




      cl_ord = sapply(1:k, function(cl) { mean(cor(eg, grps[[cl]])) - mean(cor(lg, grps[[cl]]))  })

      reord = order(order(order(cl_ord))[km$cluster])

    }
    else if (clus_method == "sort_by_cor") {
      reord = order(colMeans(cor(eg, dm_n)) - colMeans(cor(lg, dm_n)))
    }
    else if (clus_method == "hclust") {
      hc = hclust(dist(dcor), method=hc_method)
      dd = as.dendrogram(hc)
      reord.dd = reorder(dd, ord_vec, mean)
      reord.hc = as.hclust(reord.dd)
      reord = reord.hc$order
    }
  }
  d$cl_ord = reord
  
  dcor_ord = dcor[reord, reord]

  fn_pref = sprintf("%s/domain_cov_corr_karyoLT%.2f_mincov%.2f-%.2f_%s_by_cl_%s_gr%d_%s%s", odir, max_chr_cov_enr, mean_kb_cov_range[1], mean_kb_cov_range[2], switch(clus_method, kmeans=sprintf("kmeans_%dk", k), sort_by_cor="sort_by_cor", hclust=sprintf("hclust_%s", hc_method)), order_by, group_size, ifelse(ord_only, "ord_only_", ""), name)

  if (plot_fig) {
    message(sprintf("plotting using %d of %d cells...", length(use_cells), length(sch_good_cells)))
    png(paste0(fn_pref, "_cells_to_use.png"), 300, 300)
    plot(mv$ord, mv$early_f, pch=19, cex=0.5, xlab="position", ylab="repli-score")
    points(mv[use_cells, 'ord'], mv[use_cells, 'early_f'], pch=19, cex=0.6, col='red')
    dev.off()

    
    png(paste0(fn_pref, ".png"), ncol(dcor)+50, 1.5 * (ncol(dcor)+50))
    nfigs = 4 #ifelse(file.exists(ord_pc_ifn), 4, 3)
    layout(matrix(1:nfigs, nfigs, 1), 1, c(4, rep(1, times=nfigs-1)))
    par(mar=c(1,4,1,1))
    image(pmin(pmax(dcor_ord, zlim[1]), zlim[2]), zlim=zlim, col=colorRampPalette(colspec)(1000), xaxt='n', yaxt='n')
    if (group_size > 0) {
      abline(v=seq(1, ncol(dcor_ord), by=group_size)/ncol(dcor_ord), lwd=2)
    }
    
    image(1:ncol(dcor), seq(min(d$tor, na.rm=T), max(d$tor, na.rm=T), length=50), matrix(0, ncol(dcor), 50), col=NA, xaxt='n', ylab='ToR')
    points(d[as.character(colnames(dcor_ord)), 'tor'], pch=19, cex=2)
    grid(col='black', lwd=0.5)
    if (group_size > 0) {
      abline(v=seq(1, ncol(dcor_ord), by=group_size), lwd=2)
    }
    abline(h=0, col='red', lwd=4, lty=2)

    image(1:ncol(dcor), seq(min(a_s$trans_A_score, na.rm=T), max(a_s$trans_A_score, na.rm=T), length=50), matrix(0, ncol(dcor), 50), col=NA, xaxt='n', ylab='trans A-score')
    points(a_s[as.character(colnames(dcor_ord)), 'trans_A_score'], pch=19, cex=2)
    grid(col='black', lwd=0.5)
    if (group_size > 0) {
      abline(v=seq(1, ncol(dcor_ord), by=group_size), lwd=2)
    }
    abline(h=0, col='red', lwd=4, lty=2)

    #if (file.exists(ord_pc_ifn)) {
    #  v = rep(NA, ncol(dcor))
    #  names(v) = colnames(dcor_ord)
      
      #v[rownames(eord)] = eord$n
      #image(1:ncol(dcor), 1:max(v, na.rm=T), matrix(0, ncol(dcor), max(v, na.rm=T)), col=NA, xaxt='n', ylab='ExtOrd')
    ##   v[rownames(eord)] = eord$pA
    ##   image(1:ncol(dcor), seq(min(v, na.rm=T), max(v, na.rm=T), len=length(v)), matrix(0, ncol(dcor), length(v)), col=NA, xaxt='n', ylab='ExtOrd')
      
    ##   points(v, pch=19, cex=2)
    ##   grid(col='black', lwd=0.5)
    ##   if (group_size > 0) {
    ##     abline(v=seq(1, ncol(dcor_ord), by=group_size), lwd=2)
    ##   }
    ## }
    
    dc = gsub("chr", "", levels(d$chrom))
    image(1:ncol(dcor), seq_along(dc), matrix(0, ncol(dcor), 22), col=NA, xaxt='n', yaxt='n')
    axis(2, seq_along(dc), dc)
    y = as.numeric(d[colnames(dcor_ord), 'chrom'])
    points(y, pch=19, cex=2)
    grid(col='black', lwd=0.5)
    if (group_size > 0) {
      abline(v=seq(1, ncol(dcor_ord), by=group_size), lwd=2)
    }
    dev.off()
    
    plot_legend(paste0(fn_pref, "_leg.png"), zlim=zlim, crp=colorRampPalette(colspec)(1000), name="corr")
  }
  #png(sprintf("%s/domain_cov_corr_mincov%.2f_k%d_tor_boxplot.png", sch_fig_dir, min_mean_kb_cov, k), k*15, 300)
  #boxplot(tor ~ cl_ord, d, ylab='ToR')
  #dev.off()

  if (group_size > 0) {


    if (use_method == "downsample_by_g_cont") {
      dcov_ds = sch_domain_raw_marg_cov_ds[, colnames(dcor_ord)]
      missing_cells = matrix(NA, nrow(mv) - nrow(dcov_ds), ncol(dcov_ds))
      rownames(missing_cells) = rownames(mv)[!is.element(rownames(mv), rownames(dcov_ds))]
      colnames(missing_cells) = colnames(dcov_ds)
      dcov_ds = rbind(dcov_ds, missing_cells)
      dcov = dcov_ds[rownames(mv), ]
    }
    else {
      dcov_r = sch_domain_raw_marg_cov[rownames(mv), colnames(dcor_ord)]
      dcov = dcov_r / rowSums(dcov_r)
    }

    dcov_s = rollapply(t(dcov), width=group_size, by=group_size, FUN=mean, partial=T, na.rm=T)
    g1_cells = rownames(mv)[mv$group == 2]
    dcov_sn = log2(dcov_s / rowMeans(dcov_s[, g1_cells], na.rm=T))
    dcov_sns = t(rollapply(t(dcov_sn), width=sw, by=1, FUN=mean, partial=T, na.rm=T))

    nc = ceiling(sqrt(nrow(dcov_sns)))
    nr = ceiling(nrow(dcov_sns)/nc )

    if (plot_fig) { 
      png(sprintf("%s_sw%d.png", fn_pref, sw), nc*150+20, nr*120+20)
      par(mfrow=c(nr, nc))
      par(mar=c(2,2,1,1))

      apply(dcov_sns, 1, function(v) {
        plot(v, type='l', col='blue', ylim=range(dcov_sns, na.rm=T));
        apply(dcov_sns, 1, lines, col='darkgray');
        lines(v, col='blue');
        grid(col='black', lwd=0.5) })
      
      dev.off()
    }

    d$rep_bin = NA
    d[colnames(dcor_ord), 'rep_bin'] = rep(1:nrow(dcov_sns), each=group_size)[1:ncol(dcov)]
  }

  list(dcor_ord=dcor_ord, d=d, km=km, dcor=dcor, reord=reord, dm_n=dm_n)
}


###############
# 
calc_domains_ab_frac <- function(tn=pool_tn, d=NULL, min_cis=2e6, vfunc="weighted.sum", cis=F, calc_frac_a=T, rebuild=F)
{
  if (is.null(d)) {
    d = sch_get_ab_tagged_domains(raw_tn=tn, vfunc=vfunc, rebuild=rebuild)
  }

  ints = gintervals.2d.all()
  if (cis) {
    ints = filter(ints, as.character(chrom1) == as.character(chrom2))
  }
  else { 
    ints = filter(ints, as.character(chrom1) != as.character(chrom2))
  }
  d_chroms = unique(as.character(d$chrom))

  ints = ints[is.element(as.character(ints$chrom1), d_chroms) & is.element(as.character(ints$chrom2), d_chroms), ]
  
  commands = paste0(".calc_domains_ab_frac_for_chrom_pairs(d=d, tn=\"", tn, "\", chr1=\"", ints$chrom1, "\", chr2=\"", ints$chrom2, "\", min_cis=min_cis, vfunc=vfunc, calc_frac_a=calc_frac_a)", collapse=",")
  
  res = eval(parse(text=paste("gcluster.run(",commands,")")))
  stopifnot(sum(unlist(lapply(res, function(r) typeof(r$retv))) != "list") == 0)
  res_s = do.call("rbind", lapply(res, function(r) r$retv))

  if (calc_frac_a) {
    stats = as.data.frame(summarize(group_by(res_s, chrom1, start1), A=sum(A), B=sum(B)))
    stats$pA = stats$A / (stats$A + stats$B)
  }
  else {
    stats = as.data.frame(summarize(group_by(res_s, chrom1, start1), wA=sum(wA), n=sum(n)))
    stats$pA = stats$wA / stats$n
  }
  
  list(d=d, stats=stats)
}

####
.calc_domains_ab_frac_for_chrom_pairs <- function(d, tn=pool_tn, chr1="chr1", chr2="chr2", min_cis=2e6, vfunc="weighted.sum", calc_frac_a=T)
{
  stopifnot(calc_frac_a && is.element("ab", colnames(d)) || !calc_frac_a && is.element("pA", colnames(d)))

  fstr = ifelse(calc_frac_a, 'ab', 'pA')
            
  library(plyr)
  gvtrack.create("tnv", tn, vfunc)
  d1 = d[as.character(d$chrom) == chr1, ]
  d2 = d[as.character(d$chrom) == chr2, ]
  
  comb = expand.grid(1:nrow(d1), 1:nrow(d2))
  ints = gintervals.2d(d1$chrom[comb[,1]], d1$start[comb[,1]], d1$end[comb[,1]], d2$chrom[comb[,2]], d2$start[comb[,2]], d2$end[comb[,2]])

  if (chr1 == chr2) {
    x = gextract("ifelse(is.na(tnv), 0, tnv)", intervals=ints, iterator=ints, colnames='count', band=c(-1e9, -min_cis))
    x_non_trim1 = data.frame(chrom1=x$chrom1, start1=ints[x$intervalID, 'start1'], end1=ints[x$intervalID, 'end1'], chrom2=x$chrom1, start2=ints[x$intervalID, 'start2'], end2=ints[x$intervalID, 'end2'], count=x$count, intervalID=x$intervalID, d2=d2[comb[x$intervalID, 2], fstr])
    x_non_trim2 = data.frame(chrom1=x$chrom1, start1=ints[x$intervalID, 'start2'], end1=ints[x$intervalID, 'end2'], chrom2=x$chrom1, start2=ints[x$intervalID, 'start1'], end2=ints[x$intervalID, 'end1'], count=x$count, intervalID=x$intervalID, d2=d2[comb[x$intervalID, 1], fstr])

    x = rbind(x_non_trim1, x_non_trim2)

  }
  else {
    x = gextract("ifelse(is.na(tnv), 0, tnv)", intervals=ints, iterator=ints, colnames='count')
    x$d2 = d2[comb[,2], fstr]

  }

  if (calc_frac_a) {
    res = ddply(x, .(chrom1, start1), function(v) { tapply(v$count, v$d2, sum) } )
  }
  else {
    res = ddply(x, .(chrom1, start1), function(v) { data.frame(wA=sum(v$count * v$d2, na.rm=T), n=sum(v$count)) } )
  }

  res
}

################
#
calc_domains_a_score <- function(tn="scell.nextera.pool_good_hyb_2i_all_es",  vfunc="weighted.sum", d=NULL, use_trans=T, rebuild=F)
{
  ifn = sprintf("%s/%s_a_score_%s_%dDoms_%s.rdata", sch_rdata_dir, tn, vfunc, ifelse(is.null(d), 0, nrow(d)), ifelse(use_trans, "trans", "cis"))

  if (!file.exists(ifn) | rebuild) {
    message(sprintf("calc fraction of contacts with A domains using %s ...", tn))
    rd = calc_domains_ab_frac(tn=tn, d=d, cis=F, vfunc=vfunc, calc_frac_a=T, rebuild=rebuild)
    rd2 = merge(rd$d, rd$stats[, c('chrom1', 'start1', 'A', 'B', 'pA')], by.x=c('chrom', 'start'), by.y=c('chrom1', 'start1'))
    rownames(rd2) = paste(rd2$chrom, rd2$start, sep="_")

    message(sprintf("calc weigthed fraction of contacts with A domains using %s ...", tn))
    dd = calc_domains_ab_frac(tn=tn, d=rd2, cis=!use_trans, vfunc=vfunc, calc_frac_a=F, rebuild=rebuild)
    dd = dd$stats
    rd2[paste(dd$chrom1, dd$start1, sep="_"), ifelse(use_trans, 'trans_A_score', 'cis_A_score')] = dd$pA

    write.table(rd2, ifn, quote=F, sep="\t")
  }
  
  read.table(ifn, header=T)
}

############
#
sch_get_ab_tagged_domains <- function(raw_tn=pool_tn, vfunc="weighted.sum", by_ref_trans_a_score_tn=NULL, rebuild=F)
{
  ifn = sprintf("%s/%s_ab_tagged_domains_%s_by_ref_as_%s.rdata", sch_rdata_dir, raw_tn, vfunc, ifelse(is.null(by_ref_trans_a_score_tn), "TransCluster", by_ref_trans_a_score_tn))

  if (!file.exists(ifn) | rebuild) {
    ins_track = sprintf("%s_ins_%ds", raw_tn, ins_scale)
    km = cluster_pool_hic_domains_by_inter_domain_contacts(tn=raw_tn, ins_track=ins_track, vfunc=vfunc, do_plot=F)
    d = km$d
    if (mean(d[d$cluster == 1, 'tor'], na.rm=T) > mean(d[d$cluster == 2, 'tor'], na.rm=T)) { 
      ab_tags = c('A', 'B')
    }
    else {
      ab_tags = c('B', 'A')
    }
    d$ab = ab_tags[d$cluster]

    if (!is.null(by_ref_trans_a_score_tn)) {
      ref_as = calc_domains_a_score(by_ref_trans_a_score_tn, vfunc=vfunc, rebuild=rebuild)
      
      nn = gintervals.neighbors.wrap(d, ref_as, maxneighbors=100, maxdist=0, na.if.notfound=T)
      nn$start_2 = pmax(nn$start_1, nn$start_2)
      nn$end_2 = pmin(nn$end_1, nn$end_2)
      nn$len_2 = nn$end_2 - nn$start_2
      
      y = ddply(nn, .(chrom_1, start_1, len_1, ab_1), function(v) { data.frame(cov_2=sum(v$len_2), mean_as=sum(v$trans_A_score * v$len_2)/sum(v$len_2)) } )
      
      ac = hist(y[y$ab_1 == 'A', 'mean_as'], breaks=seq(0, 1, by=0.005), plot=F)
      bc = hist(y[y$ab_1 == 'B', 'mean_as'], breaks=seq(0, 1, by=0.005), plot=F)
      
      ab_cutoff = ac$mids[which.max(cumsum(bc$density) - cumsum(ac$density))]
      
      d$orig_ab = d$ab
      
      d[paste(y$chrom_1, y$start_1, sep="_"), 'ab'] = ifelse(y$mean_as >= ab_cutoff, 'A', 'B')
      
    }
    write.table(d, ifn, quote=F, sep="\t")
  }
    
  read.table(ifn, header=T)
}

########################################################
#
# Distribution of TADs A-association across single cells
#
########################################################
calc_domains_weighted_a_score_across_cells <- function(ref_tn="scell.nextera.pool_good_hyb_2i_all_es", min_cis_dist=2e6, ref_d=NULL, vfunc="weighted.sum", min_domain_contacts=4, min_unique_domains_contacted=2, use_frac_c=F, rebuild=F)
{
  if (is.null(ref_d)) {
    if (pool_tn == ref_tn) {
      ref_d = calc_domains_a_score(tn=ref_tn, vfunc=vfunc, rebuild=rebuild)
    }
    else {
      ref_d = sch_get_ab_tagged_domains(raw_tn=pool_tn, vfunc=vfunc, by_ref_trans_a_score_tn=ref_tn, rebuild=rebuild)
      ref_d = calc_domains_a_score(tn=pool_tn, vfunc=vfunc, d=ref_d, rebuild=rebuild)
      
    }
  }
  
  nms = gtrack.ls(sch_track_base)
  #nms = nms[1:10]
  
  ofn = sprintf("%s/domains_weighted_a_score_by_%s_%s_minTotal%d_minContDoms%d.Rdata", sch_rdata_dir, ref_tn, vfunc, min_domain_contacts, min_unique_domains_contacted)
  
  if (rebuild | !file.exists(ofn)) {
    commands = paste0(".sch_calc_domains_weighted_a_score(nm=\"", nms, "\", d=ref_d, min_cis_dist=min_cis_dist, min_domain_contacts=min_domain_contacts, min_unique_domains_contacted=min_unique_domains_contacted, use_frac_c=use_frac_c)", collapse=",")

    res = eval(parse(text=paste("gcluster.run(",commands,")")))
    
    mu_A = matrix(NA, nrow(ref_d), length(nms))
    rownames(mu_A) = paste(ref_d$chrom, ref_d$start, sep="_")
    colnames(mu_A) = nms
    
    var_A = total_A = n_dom_A = ds_mu_A = mu_A
    for (i in seq_along(res)) {
      rv = res[[i]]$retv
      stopifnot(typeof(rv) == 'list')

      if(length(rv) > 0) {
        row_names = paste(rv$chrom_2, rv$start_2, sep="_")
        mu_A[ row_names, nms[i]] = rv$mu_A
        ds_mu_A[ row_names, nms[i]] = rv$ds_mu_A
        var_A[ row_names, nms[i]] = rv$var_A
        total_A[ row_names, nms[i]] = rv$total
        n_dom_A[ row_names, nms[i]] = rv$n_domains
      }
    }
    save(mu_A, ds_mu_A, var_A, total_A, n_dom_A, file=ofn)
  }
  else {
    load(ofn)
  }
  
  list(mu_A=mu_A, var_A=var_A, total_A=total_A, n_dom_A=n_dom_A, ds_mu_A=ds_mu_A, ref_d=ref_d)
}

######
.sch_calc_domains_weighted_a_score <- function(nm, d, min_cis_dist=2e6, min_domain_contacts=4, min_unique_domains_contacted=2, use_frac_c=F)
{
  library(plyr)
  library(dplyr)

  res = list()

  x = gextract(nm, intervals=gintervals.2d.all(), band=c(-1e9, -min_cis_dist), colnames='count')

  
  if (!is.null(x)) {
    x1 = data.frame(chrom=x$chrom1, start=x$start1, end=x$end1)
    x2 = data.frame(chrom=x$chrom2, start=x$start2, end=x$end2)

    n1 = gintervals.neighbors.wrap(x1, d, maxdist=0, na.if.notfound=T)
    n2 = gintervals.neighbors.wrap(x2, d, maxdist=0, na.if.notfound=T)


    if (use_frac_c) {
      n1$s1 = n1$pA
      n1$s2 = n2$pA
    }
    else {
      n1$s1 = n1$trans_A_score
      n1$s2 = n2$trans_A_score
    }
    
    n1$id2 = paste0(n2$chrom_2, n2$start_2)
    
    n12 = filter(n1, !is.na(s1) & !is.na(s2))

    if (nrow(n12) > 0) {
      res = as.data.frame(summarize(group_by(n12, chrom_2, start_2), total=length(s1), n_domains=length(unique(id2)), mu_A=mean(s2), var_A=var(s2)))
      res = filter(res, total >= min_domain_contacts & n_domains >= min_unique_domains_contacted)

   
      res2 = as.data.frame(ddply(n12, .(chrom_2, start_2), function(v) { vu = unique(v[, c('id2', 's2')]); data.frame(ds_mu_A=ifelse(nrow(vu) >= 2, mean(vu[sample(1:nrow(vu), 2, F), 's2']), NA)) }))

      if (nrow(res) > 0) {
        res = merge(res, res2)
      }
      
    }
  }
  res
}

###############
#
calc_pc_pc_total_contacts <- function(tn, cis=F, pc_ifn=sprintf("%s/domain_ord_by_mean_g4_A_assoc.txt", sch_table_dir), min_cis=2e6, vfunc="weighted.sum", rebuild=F)
{
  pc = read.table(pc_ifn, header=T)
  stopifnot(is.element("pc", colnames(pc)))
  
  ifn = sprintf("%s/pc_pc_total_contacts_%s_%s_%dd_%s.rdata", sch_rdata_dir, tn, ifelse(cis, sprintf("cis_min%s", n2str(min_cis, 0)), "trans"), nrow(pc), vfunc)

  if (rebuild | !file.exists(ifn)) {
    ints = gintervals.2d.all()
    if (cis) {
      ints = filter(ints, as.character(chrom1) == as.character(chrom2))
    }
    else { 
      ints = filter(ints, as.character(chrom1) != as.character(chrom2))
    }
    pc_chroms = unique(as.character(pc$chrom))
    
    ints = ints[is.element(as.character(ints$chrom1), pc_chroms) & is.element(as.character(ints$chrom2), pc_chroms), ]
    
    commands = paste0(".calc_pc_pc_total_contacts_for_chrom_pair(tn=tn, d=pc, chr1=\"", ints$chrom1, "\", chr2=\"", ints$chrom2, "\", min_cis=min_cis, vfunc=vfunc)", collapse=",")
    
    res = eval(parse(text=paste("gcluster.run(",commands,")")))
  
    stopifnot(sum(unlist(lapply(res, function(r) typeof(r$retv))) != "double") == 0)
    stats = Reduce("+", lapply(res, function(r) r$retv))

    save(stats, file=ifn)
  }

  load(file=ifn)
  
  stats
}

###
.calc_pc_pc_total_contacts_for_chrom_pair <- function(tn, d, chr1="chr1", chr2="chr2", min_cis=2e6, vfunc="weighted.sum")
{
  library(dplyr)

  n_pc = max(d$pc)
  y = matrix(0, n_pc, n_pc)
  
  gvtrack.create("tnv", tn, vfunc)
  d1 = d[as.character(d$chrom) == chr1, ]
  d2 = d[as.character(d$chrom) == chr2, ]
  
  comb = expand.grid(1:nrow(d1), 1:nrow(d2))
  ints = gintervals.2d(d1$chrom[comb[,1]], d1$start[comb[,1]], d1$end[comb[,1]], d2$chrom[comb[,2]], d2$start[comb[,2]], d2$end[comb[,2]])

  if (chr1 == chr2) {
    x = gextract("ifelse(is.na(tnv), 0, tnv)", intervals=ints, iterator=ints, colnames='count', band=c(-1e9, -min_cis))
    x1 = data.frame(pc1=d1[comb[x$intervalID, 1], 'pc'], pc2=d2[comb[x$intervalID, 2], 'pc'], count=x$count)
    x2 = data.frame(pc1=d2[comb[x$intervalID, 2], 'pc'], pc2=d1[comb[x$intervalID, 1], 'pc'], count=x$count) 

    x = rbind(x1, x2)

  }
  else {
    x = gextract("ifelse(is.na(tnv), 0, tnv)", intervals=ints, iterator=ints, colnames='count')
    x = data.frame(pc1=d1[comb[x$intervalID, 1], 'pc'], pc2=d2[comb[x$intervalID, 2], 'pc'], count=x$count)
  }


  tot = as.data.frame(summarize(group_by(x, pc1, pc2), total=sum(count)))

  y[ cbind(tot$pc1, tot$pc2)] = tot$total
  

  y
}
