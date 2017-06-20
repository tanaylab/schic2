source(paste(Sys.getenv("PIPELINE_HOME"), "/map3c/TG3C/params.r", sep=""))
source(paste(Sys.getenv("PIPELINE_HOME"), "/map3c/TG3C/PipeStats.r", sep=""))

require("KernSmooth")

    
adj_decay_curve = function(adj_fn) {
	tab = read.table(adj_fn, header=T)

	d = ifelse(as.character(tab$chr1) == as.character(tab$chr2),
		pmax(0, log10(1+abs(tab$start1 - tab$start2)) - 3),
		-1)

	return(table(round(d*2)/2))
}

gtrack.2d.get_mat = function(track_nm, bait, horiz5, horiz3)
{
	options(gparam.type = "string")

	bait[1] = sub("chr", "", as.character(bait[1]))
	chr_id = sprintf("chr%s", bait[1])

#	message("extract from bait ", paste(bait, collapse=" "))
	bait_interv = bait
	bait_interv[2] = max(1,as.numeric(bait_interv[2])-horiz5)
	bait_interv[3] = min(as.numeric(bait_interv[3])+horiz3, 
		gintervals.all()$end[gintervals.all()$chrom==chr_id])
	
	message("forced bait")
	print(bait_interv)
	mat_interv = gintervals.2d(
			chroms1=c(as.character(bait_interv[1])),
			starts1=c(as.numeric(bait_interv[2])), 
			ends1=c(as.numeric(bait_interv[3])),
			chroms2=c(as.character(bait_interv[1])),
			starts2=c(as.numeric(bait_interv[2])), 
			ends2=c(as.numeric(bait_interv[3])))

#	message("done ini intervs\n")
	print(mat_interv)

	gvtrack.create("cisprof", track_nm, "weighted.sum")

	cont = gextract("ifelse(!is.nan(cisprof),cisprof,0)", intervals=mat_interv, iterator=track_nm, colnames=c("mol"))
	message("done extract")
	return(cont)

}

gtrack.2d.plot_mat_obs = function(track_nm, bait, horiz5, horiz3, bwx=10000,bwy=10000, logmode=0, binsize=1000) 
{
#	options(gparam.type = "string")
#	params = init_params(params_fn)
#	colspec = get_param_list("TG3C.mat_colors", params)
#	bwx = get_param("TG3C.mat_bandwidth", params)
#	bwy = get_param("TG3C.mat_bandwidth", params)

	colspec = c("white", "lightblue", "blue", "darkblue", "red",  "black", "yellow")


	mat = gtrack.2d.get_mat(track_nm, bait, horiz5, horiz3)
    if (!is.null(mat)) {
      shades = colorRampPalette(colspec)

      if(logmode == 0) {
		message("will call bkde")
		d2d = bkde2D(cbind(mat$start1, mat$start2),
          bandwidth=c(bwx,bwy), 
          gridsize=rep(binsize, 2))
		dens = d2d$fhat
		lowdens = quantile(dens[dens>0], 0.1, na.rm=T)
		highdens = quantile(dens[dens>0], 0.999,na.rm=T)
		message("dens lim ", lowdens, " and ", highdens)
		dens[dens>highdens] = highdens
                                        # ORIGINALLY
                                        # image(x=d2d$x1, y=d2d$x2, t(log(dens+lowdens)), col=shades(1000), zlim=c(log(lowdens*2), log(highdens)), useRaster=TRUE)
                                        #
		image(x=d2d$x1, y=d2d$x2, t(log(dens+lowdens)), col=shades(1000),
              zlim=c(log(lowdens*2), log(highdens)), useRaster=TRUE,
              xaxt="n", yaxt="n", xlab="", ylab="")
      } else {
		bwy = 0.25
		smoothScatter((mat$start1 + mat$start2)/2,
                      log2(abs(mat$start2-mat$start1) + 256),
                      colramp = shades, 
                      bandwidth=c(bwx, bwy))
      }
    }
  }

shuffle_contacts = function(mat, iter=10) {
	n = dim(mat)[1]
	rmat = cbind(mat[,1:3], mat[sample(n,n), 4:6])
	for(i in 2:iter) {
		message("shuf q1tart1", i)
		rmat1 = cbind(mat[,1:3], mat[sample(n,n), 4:6])
		rmat = rbind(rmat, rmat1)
	}
	return(rmat)
}

gtrack.2d.plot_mat_norm = function(track_nm, bait, horiz5, horiz3, bws, shades, mincov=30, reserve_bottom=0, maxy=(horiz3+horiz5)*0.45, ticx=2.5e+5) 
{
	m =  gtrack.2d.extract_mat_norm(c(track_nm), bait, 
					horiz5*2, horiz3*2, bws, mincov)
	baitx = as.numeric(bait[2])
	minx = baitx-horiz5
	maxx = baitx+horiz3
	image(x=m$x, y=m$y, m$mat, 
		xlim=c(minx, maxx),
		ylim=c(-reserve_bottom, maxy), 
		col=shades, useRaster=TRUE, xlab=NA,ylab=NA)

	grid_tics = seq(ticx * round((baitx-horiz5)/ticx), ticx * round((baitx+horiz3)/ticx), ticx)
	for(tic in grid_tics) {
		lines(c(minx, tic), c(tic-minx, 0), lty=2, col="black")
		lines(c(tic, maxx), c(0, maxx-tic), lty=2, col="black")
	}

}
gtrack.2d.plot_join_mat_norm = function(track_nms, bait, horiz5, horiz3, bws, shades, mincov=30) {
	
	m =  gtrack.2d.extract_mat_norm(track_nms, bait, horiz5, horiz3, bws, mincov)
	image(x=m$x, y=m$y, m$mat, col=shades, useRaster=TRUE, xlab=NA,ylab=NA)
}

gtrack.2d.cmp_mat_norm = function(track1_nm,track2_nm, bait, horiz5, horiz3, bws, shades, mincov=30, zlim=c(-3,3)) {
	m1 =  gtrack.2d.extract_mat_norm(c(track1_nm), bait, horiz5, horiz3, bws, mincov)
	m2 =  gtrack.2d.extract_mat_norm(c(track2_nm), bait, horiz5, horiz3, bws, mincov)
	m1$mat = m1$mat-mean(m1$mat, na.rm=T)
	m2$mat = m2$mat-mean(m2$mat, na.rm=T)
	message("quant m1 ", paste(quantile(m1$mat,na.rm=T), collapse=" "))
	message("quant m2 ", paste(quantile(m2$mat, na.rm=T), collapse=" "))
	image(x=m1$x, y=m1$y, m1$mat-m2$mat, col=shades, useRaster=TRUE, xlab=NA,ylab=NA, zlim=zlim)
}

gtrack.2d.extract_mat_norm = function(track_nms, bait, horiz5, horiz3, bws, mincov) {

	mat = gtrack.2d.get_mat(track_nms[1], bait, horiz5, horiz3)
	for(nm in track_nms[-1]) {
		mat =  rbind(mat, gtrack.2d.get_mat(nm, bait, horiz5, horiz3))
	}
	message("got ", dim(mat)[1], " contacts")
	v = gen_hic_mat_margnorm(mat, bws, mincov)
	return(v)
}

gen_hic_mat_margnorm = function(mat, bws, mincov=30) {

	minbw = bws[1]

	colspec = c("white", "lightblue", "blue", "darkblue", "red", "black", "black", "black", "yellow")

	shades = colorRampPalette(colspec)

	minx = min(mat$start1)

	maxbin = round(max(mat$start1-minx)/minbw)
	bins_full = table(c(1:maxbin,rep(1,maxbin),round((mat$start1-minx)/minbw)), c(1:maxbin, 1:maxbin, round((mat$start2-minx)/minbw)))

#compute marginals and 2D expected values
	marg = table(round(c(1:maxbin,(mat$start1-minx)/minbw, (mat$start2-minx)/minbw)))
	marg2d = marg %*% t(marg)

	n = dim(marg2d)[1]
	q1n = round(n/4)+1
	view_base_x = minx + q1n*minbw
	q3n = round(3*n/4)
	midn = q3n-q1n+1
	rect_ndx = seq((q1n-1)*n+q1n,(q3n-1)*n+q3n, n+1)
	rect_ndx = rep(rect_ndx, q1n) + rep(seq(0,(n-1)*(q1n-1), n-1), each =midn) 
	rbins = matrix(marg2d[rect_ndx], ncol=q1n)
	bins = matrix(bins_full[rect_ndx], ncol=q1n)
#	view_base_x=minx
#	rbins = marg2d
#	bins = bins_full

	row_cumbins = apply(bins, 1, cumsum)
	cumbins = apply(row_cumbins, 1, cumsum)
	row_rcumbins = apply(rbins, 1, cumsum)
	rcumbins = apply(row_rcumbins, 1, cumsum)

	smbins = bins
	rsmbins = rbins
	nx = dim(cumbins)[2]
	ny = dim(cumbins)[1]
	Tobs = mincov
	logasym = 1
	bws = bws[bws<(minbw*nx/4)]

	for(offset in round(bws/bws[1])) {
		of2 = offset*2
		oflog = round(log2(offset+1))
		if(logasym == 0) {
			oflog = offset
		}
		ofcol = offset+oflog
		bigbins = cumbins[(of2+1):ny, (ofcol+1):nx] +
			  cumbins[1:(ny-of2), 1:(nx-ofcol)] -
			  cumbins[(of2+1):ny, 1:(nx-ofcol)] -
			  cumbins[1:(ny-of2), (ofcol+1):nx] 
		rbigbins = rcumbins[(of2+1):ny, (ofcol+1):nx] +
			   rcumbins[1:(ny-of2), 1:(nx-ofcol)] -
			   rcumbins[(of2+1):ny, 1:(nx-ofcol)] -
			   rcumbins[1:(ny-of2), (ofcol+1):nx] 
		bigbins = cbind(matrix(nrow = ny-of2, ncol = oflog,data=0), 
			        bigbins, 
				matrix(nrow = ny-of2, ncol = offset, data=0))
		bigbins = rbind(matrix(nrow = offset, ncol = nx, data=0), 
			        bigbins, 
				matrix(nrow = offset, ncol= nx, data=0))
		rbigbins = cbind(matrix(nrow = ny-of2, ncol = oflog,data=0), 
			        rbigbins, 
				matrix(nrow = ny-of2, ncol = offset, data=0))
		rbigbins = rbind(matrix(nrow = offset, ncol = nx, data=0), 
			        rbigbins, 
				matrix(nrow = offset, ncol= nx, data=0))

		f = smbins<Tobs 
#		message("at offset ", offset, " sum ", sum(f), " comp ", sum(!f))
		smbins[f] = bigbins[f]
		rsmbins[f] = rbigbins[f]
	}
	logasym = 1
	f = smbins<(Tobs*0.6)
	smbins[f] = NA

	n = dim(smbins)[1]
	m = dim(smbins)[2]
	pad = round(bws[length(bws)]/bws[1])
	pad = 0
	return(list(x=view_base_x + bws[1]*(pad:(n-pad)), 
		    y=(0:(m-pad-1))*bws[1],
		    mat=log2(smbins/rsmbins)[pad:(n-pad),1:(m-pad)]))
	
}

gtrack.2d.gen_insu_prof = function(track_nm, scale, res, min_diag_d=1000, ignore_below=0)
{
	#extract - using an iter on  
	
    iter_1d = giterator.intervals(intervals=ALLGENOME[[1]], iterator=res)

	iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, 
			starts1=iter_1d$start, 
			ends1=iter_1d$end, 
			chroms2 = iter_1d$chrom, 
			starts2=iter_1d$start, 
			ends2=iter_1d$end)


	if(length(gvtrack.ls("obs_big")) == 1) {
		gvtrack.rm("obs_big")
	}
	if(length(gvtrack.ls("obs_ins")) == 1) {
		gvtrack.rm("obs_ins")
	}
	gvtrack.create("obs_big", track_nm, "weighted.sum")
	gvtrack.create("obs_ins", track_nm, "weighted.sum")
	gvtrack.iterator.2d("obs_big", 
		eshift1=scale, sshift1=-scale, 
		eshift2=scale, sshift2=-scale)
	gvtrack.iterator.2d("obs_ins", 
		eshift1=0, sshift1=-scale, 
		eshift2=scale, sshift2=0)

	message("will iter on ", dim(iter_2d)[1])
	ins = gextract("obs_big", "obs_ins", 
			gintervals.2d.all(),
			 iterator=iter_2d, band=c(-scale*2,-ignore_below))
    if (min_diag_d > ignore_below) {
      ins_diag = gextract("obs_big", "obs_ins", 
        gintervals.2d.all(),
        iterator=iter_2d, band=c(-min_diag_d,-ignore_below))

      ins$obs_big = ins$obs_big - ins_diag$obs_big
      ins$obs_ins = ins$obs_ins - ins_diag$obs_ins
    }

	message("will retrun ins with ", dim(ins)[1], " rows")

    return(ins)
}
gtrack.2d.gen_insu_track = function(track_nm, scale, res, min_diag_d=1000, new_track, description="", ignore_below=0)
{
	k_reg = 10
	prof = gtrack.2d.gen_insu_prof(track_nm, scale, res, min_diag_d, ignore_below=ignore_below)
	message("names ", paste(names(prof), collapse=","))
	names(prof)[1] = "chrom"
	names(prof)[2] = "start"
	names(prof)[3] = "end"

	gtrack.create_sparse(track=new_track, prof[,c(1,2,3)], log2(prof$obs_ins/(prof$obs_big+k_reg)), description=description)
}

gtrack.2d.get_insu_doms = function(insu_track, thresh, iterator=500)
{
	doms = gscreen(sprintf("is.na(%s) | %s > %f", insu_track, insu_track, thresh), iterator=iterator)
	return(doms)
}

gtrack.2d.get_insu_borders = function(insu_track, thresh, iterator=500)
{
	bords = gscreen(sprintf("!is.na(%s) & %s < %f", insu_track, insu_track, thresh), iterator=iterator)
	return(bords)
}

annot_domains_by_tss = function(doms, tss, annot_track, an_name="an")
{
	an = gextract(annot_track, doms, iterator=tss, colnames=c("an"))
	d_an = tapply(log2(1+an$an), an$intervalID, mean)
	doms$tmp = d_an[as.character(seq(1,nrow(doms),1))]
	names(doms)[ncol(doms)] = an_name
	return(doms)
}

gtrack.2d.gen_asym_prof = function(track_nm, scale, res, horiz=2000, min_diag_d=1000)
{
	#extract - using an iter on  
	
        iter_1d = giterator.intervals(intervals=ALLGENOME[[1]], iterator=res)

	iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, 
			starts1=iter_1d$start, 
			ends1=iter_1d$end, 
			chroms2 = iter_1d$chrom, 
			starts2=iter_1d$start, 
			ends2=iter_1d$end)

	gvtrack.create("obs_left", track_nm, "weighted.sum")
	gvtrack.create("obs_right", track_nm, "weighted.sum")
	gvtrack.create("obs_left_bg", track_nm, "weighted.sum")
	gvtrack.create("obs_right_bg", track_nm, "weighted.sum")
	gvtrack.iterator.2d(obs_left, 
		eshift1=0, sshift1=-scale, 
		eshift2=horiz, sshift2=-horiz)
	gvtrack.iterator.2d(obs_left_bg, 
		eshift1=0, sshift1=-scale, 
		eshift2=0, sshift2=-scale)
	gvtrack.iterator.2d(obs_right, 
		eshift2=scale, sshift2=0,
		eshift1=horiz, sshift1=-horiz)
	gvtrack.iterator.2d(obs_right_bg, 
		eshift1=scale, sshift1=0,
		eshift2=scale, sshift2=0)

	asym = gextract("obs_left", "obs_right",  "obs_left_bg", "obs_right_bg",
				gintervals.2d.all(), iterator=iter_2d, band=c(-1000000,-min_diag_d))

	return(asym)
}

add_domain_grid = function(borders, fromx, tox, col="black", lty=1)
{
	f = borders$start > fromx & borders$end < tox

	edges = round((borders$start[f] + borders$end[f])/2)
	n = sum(f)
	mids = round((edges[-1]+edges[-n])/2)
	height = round((edges[-1]-edges[-n])/2)
	message("num ", n)
	
	lines(c(c(rbind(edges[-n], mids)),edges[n]), c(c(rbind(0, height)),0), col=col, lty=lty)
}

hotspots_matrix = function(track_nms, hits1, hits2, horiz, binsize=500, diag_horiz=1e+7, norm_omit = 8, zl=c(-0.5,0.5), doms=NA, exclude_band=-1, adj_list=NA, smooth=1)
{
	#build 2d intervals by carteisna product
	is_sym = 0
	if(is.na(hits2)) {
		hits2 = hits1
		is_sym = 1
	}
	if(length(grep("strand", names(hits1))) == 0) {
		hits1$strand = 1
	}
	if(length(grep("strand", names(hits2))) == 0) {
		hits2$strand = 1
	}
	gvtrack.create("hits1_d", hits1, "distance")
	gvtrack.iterator("hits1_d", dim=1)
	gvtrack.create("hits2_d", hits2, "distance")
	gvtrack.iterator("hits2_d", dim=2)

	hspot_cart = giterator.cartesian_grid(hits1, expansion1=c(-horiz, horiz), 
						hits2, expansion2=c(-horiz,horiz))

	if(exclude_band == -1) {
		exclude_band = horiz*2
	}
	hspot1 = giterator.intervals("1", ALLGENOME, iterator=hspot_cart, band=c(exclude_band, diag_horiz))
	if(!is_sym) {
		hspot2 = giterator.intervals("1", ALLGENOME, iterator=hspot_cart, band=c(-diag_horiz, -exclude_band))
	}

	if(!is.na(doms)) {
		d11 = gintervals.annotate(
		    gintervals(hspot1$chrom1, hspot1$start1, hspot1$end1),
		    doms)
		d12 = gintervals.annotate(
		    gintervals(hspot1$chrom2, hspot1$start2, hspot1$end2),
		    doms)
		f_intra1 = (d11$doms.id == d12$doms.id)
		message("intra1 dom count ", sum(f_intra1), " out of ", dim(hspot1)[1])
		hspot1 = hspot1[f_intra1,]
		if(!is_sym) {
		  d21 = gintervals.annotate(
		    gintervals(hspot2$chrom1, hspot2$start1, hspot2$end1),
		    doms)
		  d22 = gintervals.annotate(
		    gintervals(hspot2$chrom2, hspot2$start2, hspot2$end2),
		    doms)
		  f_intra2 = (d21$doms.id == d22$doms.id)
		  hspot2 = hspot2[f_intra2,]
		}
	}

	#optionally filter intra domain by annot with dom id

	count = 1
	if(!is.na(adj_list)) {
		pairs = gextract("hits1_d", "hits2_d",  
						hspot1, iterator=adj_list)
		if(!is_sym) {
			pairs_2 = gextract("hits1_d", "hits2_d",  
						hspot2, iterator=adj_list)
			pairs = rbind(pairs, pairs_2)
		}
	} else {
		for(track_nm in track_nms) {
			pairs_1 = gextract("hits1_d", "hits2_d",  
							hspot1, iterator=track_nm)
			if(count == 1) {
				pairs = pairs_1
			} else {
				pairs = rbind(pairs, pairs_1)
			}
			if(!is_sym) {
				pairs_2 = gextract("hits1_d", "hits2_d",  
							hspot2, iterator=track_nm)
				pairs = rbind(pairs, pairs_2)
			}
			count = count + 1
		}
	}

	pairs = pairs[abs(pairs$start1-pairs$start2)>exclude_band,]

	mat = table(floor(pairs$hits1_d/binsize), floor(pairs$hits2_d/binsize))
	N = dim(pairs)[1]
	farx = (abs(pairs$hits1_d)>norm_omit*binsize)
	fary = (abs(pairs$hits2_d)>norm_omit*binsize)
	Nx = sum(fary)
	Ny = sum(farx)
	message("Total obs ", N, " max bin ", max(mat))
	matx = table(floor(pairs$hits1_d[fary]/binsize))/Nx
	maty = table(floor(pairs$hits2_d[farx]/binsize))/Ny

	mat_e = matx %*% t(maty)

	message("sun ", sum(matx))
	message("filt ", sum(mat_e*N < 2), " elements");
	mat[(mat_e * N) < 2] = NA

	shades = colorRampPalette(c("white",  "gray", "darkgray", "blue", "red", "yellow", "black"))(5000)
	mat_r = log2(mat/mean(mat,na.rm=T))-log2(mat_e/mean(mat_e,na.rm=T))
	K = dim(mat)[1]
	if(smooth) {
		mat_r_s = (mat_r[-1,-1] + mat_r[-1, -K] + mat_r[-K, -1] + mat_r[-K, -K])/4
	} else {
		mat_r_s = mat_r
	}
	message("range ", min(mat_r_s), " to ", max(mat_r_s))
	image(mat_r_s, col=shades, zlim=zl)
	return(pairs)
}

######## (yaniv) calc cis decay curve, norm as required, and optionally plot)
norm_adj_decay_curve = function(track_name, norm.by.total.obs=T, norm.by.binsize=T, min.dist=1000, plot.res=T, plot.fn=NULL, col="#00000010") {
  mat = gextract(track_name, intervals=gintervals.2d.all())
  cis = mat[as.character(mat$chrom1) == as.character(mat$chrom2) & mat$start1 <= mat$start2, ]
  d = pmax(log10(min.dist), log10(1 + abs(cis$start2 - cis$start1)))
  tab = as.data.frame(table(round(d*2)/2))
  log.dist = as.numeric(as.vector(tab$Var1))
  tab$obs = tab$Freq
  if (norm.by.total.obs) {
    tab$Freq = tab$Freq / sum(tab$Freq)
  }
  if (norm.by.binsize) {
    tab$Freq = tab$Freq / (10 ** log.dist * (sqrt(10) - 1))
  }

  if (plot.res) {
    if (!is.null(plot.fn)) {
      png(plot.fn)
      plot(as.vector(tab$Var1), log2(tab$Freq), type="l", main=track_name, xlab="distance (log10)", ylab="contact freq per bp (log2)", col=col)
      dev.off()      
    }
    else {
      lines(as.vector(tab$Var1), log2(tab$Freq), col=col)
    }
  }
  
  tab
}

#############
kmean.prepare.cis.decay=function(track.nm.re="scell.nextera.NXT_", log.dists=seq(3, 8.5, by=0.5), min.adj=10000, min.cis.frac=0.8, verbose=F)
{
  tlist = get.tracks.stats(track.nm.re)
  tlist = tlist[ tlist$n.adj >= min.adj & tlist$cis.frac >= min.cis.frac, 'track']
  if (verbose) cat(sprintf("Processing %d cells (#adj >= %d and cis-fraction >= %f)\n", length(tlist), min.adj, min.cis.frac))
  m = matrix(nrow=length(tlist), ncol=length(log.dists), data=0)
  colnames(m) = log.dists
  rownames(m) = tlist

  for (tn in tlist) {
    if (verbose) cat(paste(tn, "...\n"))
    x = norm_adj_decay_curve(tn, norm.by.total.obs=T, norm.by.binsize=F, min.dist=1000, plot.res=F)
    m[tn, x$Var1] = x$Freq
  }
  m
}

kmean.clus.cis.decay=function(track.nm.re="scell.nextera.NXT_", k=5, log.dists=seq(3, 8.5, by=0.5), min.adj=10000, min.cis.frac=0.8, cnt.q.select.range=c(0.75, 0.9), verbose=F)
{
  m  = kmean.prepare.cis.decay(track.nm.re, log.dists=log.dists, min.adj=min.adj, min.cis.frac=min.cis.frac, verbose=verbose)
  kmeans(m, centers=select.centroids(m, k, q.select.range=cnt.q.select.range), iter.max=50, trace=T)
}

############
select.centroids=function(x, k, q.select.range=c(0.75,0.9)) {
  centers = matrix(nrow=k, ncol=ncol(x), 0)
  ind = round(runif(1, min=1, max=nrow(x)))
  centers[1,] = x[ind,]
  x = x[-ind,]
  inds = c(ind)
  
  k = k - 1
  #cat(paste("chose index", ind, ", now x has", nrow(x), " items left\n"))
  while (k > 0) {
    dists = data.frame()
    for (i in seq_along(nrow(centers)-k)) {
      dists = rbind(dists, colSums((t(x) - centers[i,]) ^ 2))
    }
    dists = apply(dists, 2, function(x) min(x))
    s.d = sort(dists, index.return=T)
    ind = s.d$ix[round(runif(1, min=q.select.range[1]*nrow(x), max=q.select.range[2]*nrow(x)))]
    centers[nrow(centers)-k+1,] = x[ind,]
    
    x = x[-ind,]
    #cat(paste("chose index", ind, ", now x has", nrow(x), " items left\n"))
    k = k - 1
    inds = c(inds, ind)
  }

#  cat(paste(inds), ":")
  centers
}


#########################
track.trans.pair.contact.corr=function(track.name, chrom.pair.min.obs=30, n.min.pairs=20)
{
  vtn = paste(track.name, ".sum", sep="")
  gvtrack.create(vtn, track.name, "weighted.sum")

  trans.chrom.int = gintervals.2d.all()
  trans.chrom.int = trans.chrom.int[as.character(trans.chrom.int$chrom1) != as.character(trans.chrom.int$chrom2),]

  obs = gextract(vtn, intervals=trans.chrom.int)
  s.obs = split(obs, obs[, c("chrom1", "chrom2")])

  res = c()

  for (x in s.obs) {
    if (nrow(x) >= chrom.pair.min.obs) {
      res = c(res, cor(x$start1, x$start2))
    }
  }

  if (length(res) < n.min.pairs) {
    return(NA)
  }
  else {
    return(res)
  }
}

##########################
create.marginal.intervals=function(width, chroms=NULL)
{
  fc = gintervals.all()

  if (!is.null(chroms)) {
    if (sum(grep("chr", chroms[1])) == 0) {
      chroms = paste("chr", chroms, sep="")
    }
    fc = fc[is.element(fc$chrom, chroms),]
  }

  ints = data.frame()
  for (i in 1:nrow(fc)) {
    starts = seq(fc[i, 'start'], fc[i, 'end'], by=width)
    ends = pmin(starts + width, fc[i, 'end'])
    ints = rbind(ints, data.frame(chrom=fc[i, 'chrom'], s1=starts, e1=ends, s2=fc[i, 'start'], e2=fc[i, 'end']))
  }
  gintervals.2d(chroms1=ints$chrom, starts1=ints$s1, ends1=ints$e1, starts2=ints$s2, ends2=ints$e2)
}

########################
get.marginal=function(tn, width=10^6, chroms=c(1:19, 'X'))
{
  tns = paste(tn, "area", sep="_")
  gvtrack.create(tns, tn, "area")
  ints = create.marginal.intervals(width, chroms=chroms)
  gextract(tns, intervals=ints, iterator=ints)
}


### Pedro ###

gtrack.2d.fixedbin_marginal <- function(track, bin_size)
{
        gvtrack.create('fixbin_marg_obs', track, "weighted.sum") #double check?
        iter_1d = giterator.intervals(intervals=gintervals.all(), iterator=bin_size)
        iter_2d = gintervals.2d(iter_1d$chrom, iter_1d$start, iter_1d$end)
        binned_marginals = gextract('fixbin_marg_obs', intervals=iter_2d, iterator=iter_2d)
        binned_marginals = binned_marginals[,c(1:3,7)]
        names(binned_marginals) = c('chrom','start','end','fixbin_marg_obs')
        #TODO:optionally create a track?
        return(binned_marginals)
}


#### functions to get cisdecay from 2d tracks
get_cis <- function(tracks){
	cis_all <- gtrack.2d.cis_decay(tracks[1])
	cis_all <- cis_all[,c(2,1)]
	colnames(cis_all)[ncol(cis_all)] <- tracks[1] ; 
	for(i in tracks[-1]){
		message(i)
		c_foc = gtrack.2d.cis_decay(i)
		cis_all <- cbind(cis_all, c_foc[,1]); 
		colnames(cis_all)[ncol(cis_all)] <- i ; 
	}
	return(cis_all)	
}

gtrack.2d.cis_decay <- function(track_nm, intervals=ALLGENOME, mode="sum")
{
		#TODO : the current pipeline does not filter out non-digested pairs, therefore
		# the maps show an unbalanced signal at the close. Need to change this fn so
		# the close contacts are not taken into account
		if(mode=="sum"){
			ret_ind = 5
		}else{
			if(mode == "tot"){
				ret_ind = 1
			}else{
			stop("mode should be set to 'sum' or 'tot'")
			}	
		}
		cis <- gsummary(track_nm, intervals=intervals, band=c(2^0, 2^10))[ret_ind]
		
		for(i in 10:28){
			cis <- rbind(cis, gsummary(track_nm, intervals=intervals, band=c(2^i, 2^(i+1)))[ret_ind])
		}
		cis[is.na(cis)] <- 0
		cis <- cbind(cis, 2^(10:29))
		colnames(cis) <- c("sum", "bin")
		rownames(cis) <- NULL
		cis[,1] <- cis[,1]/sum(cis[,1], na.rm=T)
		return(as.data.frame(cis))
}



# generate Ren's directionality index
gtrack.2d.gen_DI_prof = function(track_nm, scale=2e6, res=40e3, min_diag_d=1000)
{
  iter_1d = giterator.intervals(intervals=ALLGENOME[[1]], iterator=res)

  iter_2d = gintervals.2d(chroms1 = iter_1d$chrom, 
    starts1=iter_1d$start, 
    ends1=iter_1d$end, 
    chroms2 = iter_1d$chrom, 
    starts2=iter_1d$start, 
    ends2=iter_1d$end)

  browser()
  if(length(gvtrack.ls("obs_up")) == 1) {
    gvtrack.rm("obs_up")
  }
  if(length(gvtrack.ls("obs_down")) == 1) {
    gvtrack.rm("obs_down")
  }

  gvtrack.create("obs_up", track_nm, "weighted.sum")
  gvtrack.create("obs_down", track_nm, "weighted.sum")
  
  gvtrack.iterator.2d("obs_up", sshift1=-scale, eshift1=-res)
  gvtrack.iterator.2d("obs_down", sshift2=res, eshift2=scale)
  
  message("will iter on ", dim(iter_2d)[1])
  ins = gextract("obs_up", "obs_down", gintervals.2d.all(), iterator=iter_2d, band=c(-scale-res,-min_diag_d))
  
  message("will retrun ins with ", dim(ins)[1], " rows")

  return(ins)
}

