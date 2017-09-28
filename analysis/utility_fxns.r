#########
fig_png <- function(filename, width=480, height=480)
{
  if (plot_type == "png") {
    png(filename, width=width*fig_factor, height=height * fig_factor, pointsize=pointsize, res=fres)
  }
  else if (plot_type == "ps") {
    #postscript(file=sub("png$", "eps", filename), width=(width*fig_factor) / fres, height=(height * fig_factor) / fres, horizontal = FALSE, onefile = FALSE, paper = "special", pointsize=pointsize)
    cairo_ps(file=sub("png$", "eps", filename), width=(width*fig_factor) / fres, height=(height * fig_factor) / fres, onefile = FALSE, pointsize=pointsize)

  }
}

#########
image_table <- function(x, title=character(0), col=c(), zero_col=NULL, colrange=c(), show_xlab=T, show_ylab=T, xlabels=c(), ylabels=c(), xlabels_at=c(), ylabels_at=c(), grid_at=c(), grid_col="black", cex=1, nshades=1000, square=F, colorkey=T, cor=F)
{
    if (show_xlab) {
        if (is.null(xlabels)) {
            xlabels <- colnames(x)
        }

        if (is.null(xlabels_at)) {
            xlabels_at <- c(1:length(xlabels))
        }
    }
    if (show_ylab) {
        if (is.null(ylabels)) {
            ylabels <- rownames(x)[nrow(x):1]
        }

        if (is.null(ylabels_at)) {
            ylabels_at <- c(1:length(ylabels))
        }
    }

    # color stuff
    # sequential favorites
    pubu <- c("#fff7fb","#ece7f2","#d0d1e6","#a6bddb","#74a9cf","#3690c0","#0570b0","#045a8d","#023858")
    pubugn <- c("#fff7fb","#ece2f0","#d0d1e6","#a6bddb","#67a9cf","#3690c0","#02818a","#016c59","#014636")
    bupu <- c("#f7fcfd","#e0ecf4","#bfd3e6","#9ebcda","#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b")
    gnbu <- c("#f7fcf0","#e0f3db","#ccebc5","#a8ddb5","#7bccc4","#4eb3d3","#2b8cbe","#0868ac","#084081")

    ylgnbu <- c("#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58")
    reds <- c("#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c","#cb181d","#a50f15","#67000d")
    oranges <- c("#fff5eb","#fee6ce","#fdd0a2","#fdae6b","#fd8d3c","#f16913","#d94801","#a63603","#7f2704")
    ylorbr <- c("#ffffe5","#fff7bc","#fee391","#fec44f","#fe9929","#ec7014","#cc4c02","#993404","#662506")

    # diverging favorites
    rdbu <- c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac")
    rdylgn <- c("#a50026","#d73027","#f46d43","#fdae61","#fee08b","#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850","#006837")

    # set the color ramp, and reverse if needed
    if (is.null(col)) {
        if (!is.null(zero_col)) {
            pubu <- pubu[2:length(pubu)]
        }
        color.ramp <- colorRampPalette(pubu)(nshades)
    }
    else {
        colors = get(sub("\\*","",col))
        if (!is.null(zero_col)) {
            colors <- colors[2:length(colors)]
        }
        color.ramp <- colorRampPalette(colors)(nshades)
        if (grepl("\\*", col)) {
            color.ramp <- rev(color.ramp)
        }
    }

    if (!is.null(zero_col)) {
        color.ramp[1] <- zero_col
    }

    # if no color range specified, choose according to correlation or just min/max
    if (is.null(colrange)) {
        if (cor==T) {
            bound <- max(abs(x[!is.na(x) & x!=1]))
            color.levels <- seq(-bound, bound, length=nshades)
            colrange <- c(-bound, bound)
        }
        else {
            color.levels <- seq(min(x[!is.na(x)]), max(x[!is.na(x)]), length=nshades)
            colrange <- c(color.levels[1], color.levels[length(color.levels)])
        }
    }
    else {
        color.levels <- seq(colrange[1], colrange[2], length=nshades)
    }

    # plot and colorkey dimensions
    if (colorkey) {
        if (square) {
            map.rel.width <- 8
            pty.hm="m"
        }
        else {
            map.rel.width <- 18
            pty.hm="m"
        }
        pty.ck="m"
        layout(matrix(data=c(rep(1,map.rel.width),rep(2,2)), nrow=1, ncol=map.rel.width+2), heights=c(1,1))
    }
    else {
        pty.hm="m"
    }

    if (length(title)==0) {
        mar.heatmap <- c(2,5*cex,4*cex,0)
        mar.colorkey <- c(2,1,4*cex,4*cex)
    }
    else {
        mar.heatmap <- c(2,5*cex,10*cex,0)
        mar.colorkey <- c(2,1,10*cex,4*cex)
    }

    # matrix heatmap
    par(mar=mar.heatmap, pty=pty.hm)
    if (is.null(colrange)) {
        image(x=1:ncol(x), y=1:nrow(x),
            # z=t(x)[nrow(x):1,nrow(x):1],
            z=t(x),
            ylim=c(nrow(x),1),
            xlab="", ylab="", axes=F,
            col=color.ramp)
    }
    else {
        image(x=1:ncol(x), y=1:nrow(x),
            # z=t(x)[nrow(x):1,nrow(x):1],
            z=t(x),
            ylim=c(nrow(x),1),
            xlab="", ylab="", axes=F,
            col=color.ramp, zlim=colrange)
    }
    if (!is.null(grid_at)) {
        abline(h=grid_at, v=grid_at, col=grid_col)
    }
    box()
    if (show_xlab) {
        axis(side=3, at=xlabels_at, labels=xlabels, cex.axis=2*cex)
    }
    if (show_ylab) {
        axis(side=2, at=ylabels_at, labels=ylabels, las=HORIZONTAL<-1, cex.axis=2*cex)
    }

    if (length(title)!=0) {
        title(title, cex.main=cex*2.5)
    }

    # color key
    if (colorkey) {
        par(mar=mar.colorkey, pty=pty.ck)
        # par(pty=pty.ck)
        image(1, color.levels,
            matrix(data=color.levels, ncol=length(color.levels), nrow=1),
            col=color.ramp, xlab="",ylab="", axes=F)
        box()
        axis(4, cex.axis=cex*2, las=2)
    }
}

#########
# from stack overflow

simple_cap <- function(x)
{
    # s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(x,1,1)), substring(x,2), sep="", collapse=" ")
}

#########

# Wrapper to AT K-means clustering (impl by C) that supports NA values
TGLKMeans_wrapper <- function(data, fn, k)
{
    write.table(data, fn, sep="\t", quote=F, col.names=NA)
    system(sprintf("./bin/TGLKMeans_static %s %s euclid -allow_nas=1 &> %s.log", fn, k, fn))
    km = list()
    m.k = read.table(paste(fn, "kclust", sep="."), header=T)
    m.k$clust = m.k$clust + 1
    m.c = read.delim(paste(fn, "center", sep="."), header=F)
    km$size = tapply(m.k$id, m.k$clust, length)
    km$cluster = m.k$clust
    names(km$cluster) = m.k$id
    km$centers = m.c[, 2:ncol(m.c)]
    km
}

###################
# Input: list of cells, chromosome name, quantile to filter counts on, and which contacts to count (options: self, X2M, far and trans)
# Output: List of filtered cells and the minimum number of contacts they have on the chromosome
# Note that the contacts are doubled, so x-y and y-x are considered 2 contacts
filter_by_contacts_per_chrom <<- function(cells, chrom, quant, filter_by=c('X2M', 'far'))
{
    chrom = ifelse(length(grep("^chr", chrom)) > 0, chrom, paste("chr", chrom, sep=""))
    cols = paste(filter_by, chrom, sep="_")
    x = rowSums(sch_chrom_stat[cells, cols])
    n = quantile(x, quant) - 1
    x = x[x >= n]
    list(cells=names(x), min_contacts=n)
}


###################
vals_to_cols <- function(vals, breaks, ncols=256)
{
    min = breaks[1]
    max = breaks[length(breaks)]
    n = length(breaks)-1
    cols = rep(-1, length(vals))
    for (i in 1:n)
    {
        ind = !is.na(vals) & (breaks[i] <= vals) & (vals <= breaks[i+1])
        if (!any(ind))
            next
        # normalize to [0,1]
        cols[ind] = (vals[ind] - breaks[i]) / (breaks[i+1] - breaks[i])
        # normalize to [i*ncols,i*(ncols+1)]
        cols[ind] = (i-1)*ncols + cols[ind]*(ncols-1) + 1
        # round
        cols[ind] = round(cols[ind])
    }
    cols[cols == -1] = NA
    return (cols)
}

########################
plot_legend <- function(ofn, zlim, crp, name)
{
    fig_png(ofn, width=3*length(crp), height=150)
    image(as.matrix(seq(zlim[1], zlim[2], length=length(crp))), zlim=zlim, col=crp, xaxt='n', yaxt='n', main=name)
    axis(1, at=c(0,1), labels=zlim)
    dev.off()
}

########################
n2str <- function(n, digits=2)
{
    f    = ifelse(n >= 1e9, 1e9, ifelse(n >= 1e6, 1e6, ifelse(n >= 1e3, 1e3, 1)))
    suff = ifelse(n >= 1e9, "G", ifelse(n >= 1e6, "M", ifelse(n >= 1e3, "K", "")))

    ifelse(digits > 0, sprintf(paste("%.", digits, "f%s", sep=""), n/f, suff), paste(n/f, suff, sep=""))
}

########################## addalpha() from https://github.com/mylesmharrison/colorRampPaletteAlpha/blob/master/colorRampPaletteAlpha.R
addalpha <- function(colors, alpha=1.0) {
    r <- col2rgb(colors, alpha=T)
    # Apply alpha
    r[4,] <- alpha*255
    r <- r/255.0
    return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

############

load_table_if_exists <- function(fn, header=T, sep="\t") {
    out = NULL
    if (file.exists(fn)) {
        out = read.table(fn, header=header, sep=sep, check.names=F)
    }
    out
}


############
cyclic_rollapply <- function(x, width, by=1, func_str="mean", na.rm=T) {
    stopifnot(width %% 2 == 1)
    hsw = (width-1)/2
    if (is.vector(x)) {
        l = length(x)
        y = x[c((l - hsw + 1):l, 1:l, 1:hsw)]
    }
    else {
        l = ncol(x)
        y = t(x[,c((l - hsw + 1):l, 1:l, 1:hsw)])
    }

    if (func_str == "mean") {
        y = rollmean(y, k=width, fill=NA)
    }
    else if (func_str == "sum") {
        y = rollsum(y, k=width, fill=NA)
    }
    else if (func_str == "sd") {
        y = rollapply(y, width=width, FUN=sd, partial=T, na.rm=T)
    }
    else if (func_str == "length") {
        y = rollapply(y, width=width, FUN=length, partial=T)
    }

    if (is.vector(x)) {
        y = y[seq(hsw + 1, hsw + l,  by=by)]
    }
    else {
        y = t(y[ seq(hsw + 1, hsw + l,  by=by), ])
    }
    y
}

cyclic_rollmean <- function(x, width, by=1) {
    cyclic_rollapply(x, width, by, "mean")
}

cyclic_rollsum <- function(x, width, by=1) {
    cyclic_rollapply(x, width, by, "sum")
}

##############
# intervals manipulation (from Netta)
intervals.centers <- function(inv){
    inv[,2:3]<-floor((inv[,2]+inv[,3])/2)
    inv[,3]<-inv[,3]+1
    return(inv)
}

intervals.2d.centers <- function(inv){
  inv[,2:3]<-floor((inv[,2]+inv[,3])/2)
  inv[,3]<-inv[,3]+1
  inv[,5:6] = floor((inv[,5]+inv[,6])/2)
  inv[,6] = inv[,5] + 1
  return(inv)
}

intervals.expand <- function(inv,expansion=100){
  inv[,2]<-inv[,2]-expansion
  inv[,3]<-inv[,3]+expansion
  return(gintervals.force_range(inv))
}

intervals.2d.expand <- function(inv,expansion1, expansion2){
  inv[,2]<-inv[,2]-expansion1
  inv[,3]<-inv[,3]+expansion1
  inv[,5]<-inv[,5]-expansion2
  inv[,6]<-inv[,6]+expansion2
  return(gintervals.force_range(inv))
}


##############
# http://stackoverflow.com/questions/3318333/split-a-vector-into-chunks-in-r
chunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

##############

diverge_color <- function(colors, break_vals, nshades)

{
  ramps = list()
  ranges = c()
  for (i in 1:(length(colors)-1)) {
    ramps[[i]] = colorRampPalette(c(colors[i],colors[i+1]))
    ranges[i] = break_vals[i+1] - break_vals[i]
  }

  nshade_list = round(nshades*(ranges/(max(break_vals)-min(break_vals))))

  my_colors = c()
  for (i in 1:(length(colors)-1)) {
    my_colors = c(my_colors, ramps[[i]](nshade_list[i]))
  }
  return(my_colors)
}

###########

.plot_points_and_trend <- function(x, ofn, cols, width_per_panel=360, height_per_panel=110, sw=101, col_by_cond=T, collapse_q=c(0.005, 0.995), cex=0.4, cex.lab=1, glob_trend=NULL, glob_col=NA, glob_ylim=NULL, lcol="blue", pcol="black", vertical=T, plot_trends=F, pointsize=7, fres=120, trend_cols=mini_qualitative_colors, group_size_for_boxplot=NA, boxplot_fill='blue', show_pvals=F, boxplot_fdr=0.05, selected_cols=NULL, boxplot_comp_inds=NULL, show_pvals_as_nums=F, lwd=2, boxplot_grid=F, show_yaxis=F, outlier_col='red', inline=F, labels=NULL, at=NULL)
{

  if (!is.null(selected_cols)) {
    x = x[, selected_cols]
  }

  np = ncol(x) + ifelse(plot_trends, 1, 0) + ifelse(is.null(labels), 0, 1)

  ylims = matrix(0, np, 2)

  if (vertical) {
    if (!inline) {
      fig_png(ofn, width_per_panel, height_per_panel * np)
    }
    par(mar=c(1, 4, 0.2, 1))
    layout(matrix(1:np, np, 1))
  }
  else {
    if (!inline) {
      fig_png(ofn, width_per_panel * np, height_per_panel)
    }
    par(mar=c(4, 2, 0.2, 0.2))
    layout(matrix(1:np, 1, np))
  }

  trends = NULL
  for (i in 1:ncol(x)) {
    v = x[,i]

    if (col_by_cond) {
      col = cols
    }
    else {
      col = rep(pcol, length(v))
    }

    if (is.null(glob_ylim)) {
      qs = quantile(v, collapse_q, na.rm=T)
    }
    else {
      qs = glob_ylim
    }

    ind = v < qs[1] | v > qs[2]
    v[v < qs[1]] = qs[1]
    v[v > qs[2]] = qs[2]

    if (is.null(glob_ylim)) {
      ylim = range(v, na.rm=T)
    }
    else {
      ylim = glob_ylim
    }
    ylims[i,] = ylim

    col[ind] = outlier_col

    plot(v, pch=19, cex=cex, cex.lab=cex.lab, cex.axis=cex.lab, main="", ylab="", yaxt='n', xaxt='n', xlab="", col=col, ylim=ylim)
    if (vertical) {
      axis(2)
      title(ylab=colnames(x)[i])
    }
    else {
      if (i == 1 | show_yaxis) {
        axis(2)
      }
      title(xlab=colnames(x)[i])
    }

    if (!is.null(glob_trend)) {
      lines(glob_trend, col=glob_col, lwd=lwd)
    }
    vs = cyclic_rollmean(v, sw)
    lines(vs, col=lcol, lwd=lwd)
    grid(col='black', nx=NA, ny=NULL)

    trends = rbind(trends, vs)
  }


  if (plot_trends) {
    plot(trends[1,], type='l', ylim=range(trends), col=NA, main="", xlab="", ylab="")
    apply(cbind(trend_cols[1:nrow(trends)], trends), 1, function(y) { lines(y[-1], col=y[1], lwd=lwd) })
    grid(col='black')
    legend("topleft", legend=colnames(x), col=trend_cols[1:nrow(trends)], lwd=2, bty='n', cex=0.8, ncol=2)
  }

  if (!is.null(labels)) {
    axis(ifelse(vertical, 1, 4), labels=labels, at=at, cex.lab=cex.lab, tick=F)
  }
  
  if (!inline) {
    dev.off()
  }

  # plot data in boxplots if required
  if (!is.na(group_size_for_boxplot) ) {
    np = ncol(x)
    bp_ofn = gsub(".png", sprintf("_grp%d%s.png", group_size_for_boxplot, ifelse(show_pvals, "_pv", "")), ofn)
    if (vertical) {
      if (!inline) {
        fig_png(bp_ofn, width_per_panel, height_per_panel * np)
      }
      par(mar=c(1, 4, 0.2, 1))
      layout(matrix(1:np, np, 1))
    }
    else {
      if (!inline) {
        fig_png(bp_ofn, width_per_panel * np, height_per_panel)
      }
      par(mar=c(4, 2, 0.2, 0.2))
      layout(matrix(1:np, 1, np))
    }
    bins = ceiling( (1:nrow(x)) / group_size_for_boxplot)

    for (i in 1:ncol(x)) {
      v = x[,i]

      ylim= ylims[i,]

      v = pmin(pmax(v, ylim[1]), ylim[2])

      if (show_pvals) {
        ps = ylim[2]* 1.1
        ylim[2] = ps * 1.1
        pe = ylim[2]

        if (is.null(boxplot_comp_inds)) {
          inds = cbind(1:max(bins), (1:max(bins))+1)
        }
        else {
          inds = boxplot_comp_inds
        }

        pvals = c()
        for (r in 1:nrow(inds)) {
          pvals = c(pvals, ks.test(v[bins == inds[r,1]], v[bins == ifelse(inds[r,2] == max(bins)+1, 1, inds[r,2])])$p.value)
        }
        pvals = p.adjust(pvals, 'fdr')
        print(cbind(inds, pvals))
        if (show_pvals_as_nums) {
          pvals_lab = sapply(pvals, function(p) { sprintf("%.3e", p) })
        }
        else {
          pvals_lab = ifelse(pvals < 1e-3, '***', ifelse(pvals < 1e-2, '**', ifelse(pvals < 5e-2, '*', '-')))
        }
      }

      boxplot(split(v, bins), notch=T, ylim=ylim, xaxt='n', col=boxplot_fill, pch=19, cex=0.5)
      outliers = v == ylim[1] | v == ylim[2]
      points(bins[outliers], v[outliers], pch=19, cex=0.5, col='red')

      if (boxplot_grid) {
        grid(nx=NA, ny=NULL,  col='black', lwd=0.5)
      }

      if (vertical) {
        title(ylab=colnames(x)[i])
      }
      else {
        title(xlab=colnames(x)[i])
      }

      if (show_pvals) {
        if (is.null(boxplot_comp_inds)) {
          inds = inds[pvals <= boxplot_fdr, ]
          pvals_lab = pvals_lab[pvals <= boxplot_fdr]
        }

        bh = (pe - ps) / 4

        maxh = tapply(v, bins, max)
        maxph = apply(inds, 1, function(b) { max(maxh[seq(b[1], ifelse(b[2] == max(bins)+1, 1, b[2]))])})
        maxph = maxph + bh

        if (nrow(inds) > 0) {
          segments(inds[,1], maxh[inds[,1]]+bh, inds[,1], maxph+bh, col='red', lwd=0.5)
          segments(inds[,2], maxh[inds[,2]]+bh, inds[,2], maxph+bh, co='red', lwd=0.5)
          segments(inds[,1], maxph+bh, inds[,2], maxph+bh, col='red', lwd=0.5)
          text(rowMeans(inds), maxph+bh, pos=3, labels=pvals_lab, cex=0.8, offset=0.1, col='red', srt=30)
        }
      }
    }
    if (!inline) {
      dev.off()
    }
  }
}

.get_well_covered_scope <- function (iter=20e3, otn=pool_tn, min_cis_marg_cov_per_kb=50)
{
  i1 = giterator.intervals(intervals=gintervals.all(), iterator=iter)
  i2 = gintervals.2d(i1$chrom, i1$start, i1$end, i1$chrom, 0, -1)
  gvtrack.create("otn_area", otn, "area")

  iscope = gscreen(sprintf("otn_area >= %d", (min_cis_marg_cov_per_kb * iter) / 1e3),  intervals=i2, iterator=i2)
  scope1 = gintervals(iscope$chrom1, iscope$start1, iscope$end1)
  gintervals.union(scope1, scope1)
}

gintervals.neighbors.wrap <- function(intervals1, intervals2, maxneighbors = 1,
                          mindist = -1e+09, maxdist = 1e+09,
                          mindist1 = -1e+09, maxdist1 = 1e+09,
                          mindist2 = -1e+09, maxdist2 = 1e+09,
                          na.if.notfound = FALSE)

{
  nn = gintervals.neighbors(intervals1, intervals2, maxneighbors = maxneighbors,
    mindist = mindist, maxdist = maxdist,
    mindist1 = mindist1, maxdist1 = maxdist1,
    mindist2 = mindist2, maxdist2 = maxdist2,
    na.if.notfound = na.if.notfound)

  if (!is.null(nn) & length(intersect(colnames(intervals1), colnames(intervals2))) > 0) {
    renames = c(paste0(colnames(intervals1), '_1'), paste0(colnames(intervals2), '_2'))
    colnames(nn)[seq_along(renames)] = renames

  }

  nn
}

##############
get_cond_cols <- function(haploid=F) {
  
  cols = rep('darkgray', nrow(sch_decay_metrics))
  if (haploid) {
    cind = sch_decay_metrics$cond %in% names(sch_cond_colors)
  }
  else {
    cind = sch_decay_metrics$cond %in% grep("1CD_", names(sch_cond_colors), value=T)    
  }
  cols[cind] = unlist(sch_cond_colors[sch_decay_metrics$cond[cind]])

  cols
}

###############
partition_contacts <- function(nm, n_partitions, chroms, seed=NULL)
{
  chroms = paste0("chr", chroms)
  scope = gintervals.2d.all()
  scope = scope[is.element(scope$chrom1, chroms) & is.element(scope$chrom2, chroms), ]
  
  x = gextract(nm, intervals=scope, colnames='count')
  
  set.seed(ifelse(is.null(seed), n_partitions, seed)) 
  x$part = sample(1:n_partitions, nrow(x), replace=T)
  
  x
}
##############
# https://en.wikipedia.org/wiki/Jensen%E2%80%93Shannon_divergence#Definition

dKL <- function(p, q) {
    # TO DO if length(p) != length(q), give error..
    sum(sapply(1:length(p), function(i) {
        ifelse(p[i]==0, 0, p[i] * log(p[i]/q[i]))
    }))
}

dJS <- function(p, q) {
    m = 0.5 * (p + q)
    # 0.5 * (sum(p * log(p / m)) + sum(q * log(q / m)))
    0.5 * (dKL(p, m) + dKL(q, m))
}

dJS_dist <- function(mat) {
    apply(mat, 1, function(x) {
        apply(mat, 1, function(y) {
            foo = sqrt(dJS(x,y))
            if (is.na(foo)) { message("ah!") }
            return(foo)
            })
        })
}


######################
.update_shaman_options <- function()
{
  options(shaman.score_pal_0_col=diver_colspec[6])
  options(shaman.score_pal_pos_col=diver_colspec[7:11])
  options(shaman.score_pal_neg_col=diver_colspec[5:1])
  options(shaman.score_pal_pos_breaks=c(15, 25, 45, 60, 85))
  options(shaman.score_pal_neg_breaks=c(15, 25, 45, 60, 85))
}
