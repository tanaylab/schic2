#library(MASS)
#library(zoo)

##################################
#
# Assumes prior calls for:
#    - sch_load_tables
#
# insulation schematically: #contacts in A /  (A+B1+B2)
# store results in tables
#
#           --------
#           |   |B2/
#           | A | /
#           |___|/
#           |B1 /
#           |  /
#           | /
#           |/
#      
##################################

###########################################################################
# downsample K contacts independently in each triangle of a locus (A+B1+B2)
sch_calc_insu_ds_trian <<- function(chrom, n_ds=30, coords=NULL, coords_nm="", discard_below=1000, return_data=T, rebuild=F, scale=insu_scale, res=insu_res)
{
  sch_chrom_stat <<- read.table(sprintf("%s/chrom_stat.txt", sch_table_dir),header=T)
  
  
  system(sprintf("mkdir -p %s/insu", sch_table_dir))

  tot = list(a = NULL, b1 = NULL, b2 = NULL, ds_a = NULL, ds_b1 = NULL, ds_b2 = NULL)
  fns = list()
  for (tn in names(tot)) {
    if (length(grep("ds", tn)) > 0) {
      fns[tn] = sprintf("%s/insu/%s_insu_%d_every_%d_ds%d_%s%s.txt", sch_table_dir, chrom, scale, res, n_ds, tn, coords_nm)
    }
    else {  
      fns[tn] = sprintf("%s/insu/%s_insu_%d_every_%d_%s%s.txt", sch_table_dir, chrom, scale, res, tn, coords_nm) 
    }
  }

  if (sum(!sapply(fns, file.exists)) == 0 && !rebuild) {
    for (tn in names(tot)) {
      #message(sprintf("Reading %s", fns[[tn]]))
      tot[[tn]] = read.table(fns[[tn]], header=T)
    }
  }
  else {
    chr_coords = gintervals(chrom, 0, -1)
    scope = gintervals.2d(chrom, 0, -1, chrom, 0, -1)

    if (is.null(coords)) {
      coords = seq(chr_coords$start + scale, chr_coords$end - scale - 1, by=res)
    }
    coords = pmin(coords, chr_coords$end - 1)
	iter_2d = gintervals.2d(chroms1 = chrom, 
			starts1=coords, 
			ends1=coords + 1,
			chroms2 = chrom, 
			starts2=coords,
			ends2=coords+1)
    
    for (nm in rownames(sch_chrom_stat)) {
      message("Processing ", nm, "...")
      
      if(length(gvtrack.ls("obs_a")) == 1) {
		gvtrack.rm("obs_a")
      }
      if(length(gvtrack.ls("obs_b1")) == 1) {
		gvtrack.rm("obs_b1")
      }
      if(length(gvtrack.ls("obs_b2")) == 1) {
		gvtrack.rm("obs_b2")
      }
      gvtrack.create("obs_a", nm, "area")
      gvtrack.create("obs_b1", nm, "area")
      gvtrack.create("obs_b2", nm, "area")

      gvtrack.iterator.2d("obs_a", 
                          sshift1=-scale, eshift1=0, 
                          sshift2=0, eshift2=scale)

      gvtrack.iterator.2d("obs_b1", 
                          eshift1=-res, sshift1=-scale, 
                          eshift2=-res, sshift2=-scale)

      gvtrack.iterator.2d("obs_b2", 
                          eshift1=scale, sshift1=res, 
                          eshift2=scale, sshift2=res)

      ins = gextract("obs_a", "obs_b1", "obs_b2", scope, iterator=iter_2d, band=c(-scale*2,1))
      diag = gextract("obs_a", "obs_b1", "obs_b2", scope, iterator=iter_2d, band=c(-discard_below,1))
      st = data.frame(a = ins$obs_a - diag$obs_a, 
        b1 = ins$obs_b1 - diag$obs_b1,
        b2 = ins$obs_b2 - diag$obs_b2)

      ds_st = apply(st, 1, function(x) { if(!anyNA(x) && sum(x, na.rm=T) >= n_ds) { table(c(names(st), sample(rep(names(st), x), n_ds))) - 1 } else { rep(NA, 3) } })
      ds_st = as.data.frame(t(ds_st))
      colnames(ds_st) = names(st)
      
      tot$a = rbind(tot$a, st$a)
      tot$b1 = rbind(tot$b1, st$b1)
      tot$b2 = rbind(tot$b2, st$b2)
      
      tot$ds_a = rbind(tot$ds_a, ds_st$a)
      tot$ds_b1 = rbind(tot$ds_b1, ds_st$b1)
      tot$ds_b2 = rbind(tot$ds_b2, ds_st$b2)

    }
    for (tn in names(tot)) {
      rownames(tot[[tn]]) = rownames(sch_chrom_stat)
      colnames(tot[[tn]]) = coords    
      write.table(tot[[tn]], fns[[tn]], quote=F, sep="\t")
    }
  }
  
  if (return_data) {
    return(tot)
  }
}






#####
###get ds data a per cell
get_ds_insu_stats <- function(chroms=paste0("chr", sch_chroms), cells=sch_good_cells, n_ds=30, scale=insu_scale, res=insu_res, coords_nm="", rebuild=F) {
  m = list(a=NULL, b1=NULL, b1=NULL, ds_a=NULL, ds_b1=NULL, ds_b2=NULL)
  for (chr in chroms) {
    message("Processing chrom ", chr)
    cr = sch_calc_insu_ds_trian(chr, n_ds, coords_nm=coords_nm, return_data=T, scale=scale, res=res, rebuild=rebuild)
    for (s in names(cr)) {
      colnames(cr[[s]]) = paste(chr, colnames(cr[[s]]), sep="_")
      m[[s]] = rbind(m[[s]], t(cr[[s]][cells,]))
    }
  }
  m
}


####
get_tad_borders_per_group <- function(g_nms = c('post_m', 'g1', 'early_s', 'mid_s_g2', 'pre_m'), max_bord_dist=20e3, g_tns=NULL, scope=NULL)
{
  if (is.null(g_tns)) {
    g_tns=paste(pool_tn, "group", 1:5, g_nms, sep="_")
  }
  if (is.null(scope)) {
    scope = .get_well_covered_scope(5e3, pool_tn, 50)
  }
  borders = NULL
  bl = list()
  gbl = list()
  for (i in seq_along(g_tns)) {
    nm = g_nms[i]
    message(nm)
    bl[[nm]] = .sch_get_pool_tad_borders(input_track=g_tns[i], ins_track=sprintf("%s_ins_%ds", g_tns[i], ins_scale), scale=ins_scale, min_diag=sch_remove_near_dist)

    gbl[[nm]] = gintervals.intersect(bl[[nm]], scope)
    
    message(nrow(bl[[nm]]), " borders, out of them ", nrow(gbl[[nm]]), " well covered")
    borders = rbind(borders, bl[[nm]])
  }
  

  
  bn = gintervals.neighbors(borders, borders, maxneighbors=100, mindist=1000, maxdist=max_bord_dist, na.if.notfound=T)

  buc = gintervals.canonic(rbind(bn[is.na(bn$dist), 1:3], data.frame(chrom=bn[!is.na(bn$dist), 'chrom'], start=apply(bn[!is.na(bn$dist), grep("start", colnames(bn))], 1, min), end=apply(bn[!is.na(bn$dist), grep("start", colnames(bn))], 1, max))))
  buc = buc[ buc$end - buc$start <= max_bord_dist, ]
  buc = intervals.centers(buc)

  gbuc = gintervals.intersect(buc, scope)

  g_ins = paste0(g_tns, "_ins_", format(ins_scale, scientific=F), "s")
  b_ins_vals = gextract(g_ins, intervals=gintervals(sort(unique(gbuc$chrom))), iterator=gbuc, colnames=g_nms)
  
  list(scope=scope, buc=buc, bord=bl, gbuc=gbuc, gbord=gbl, b_ins_vals=b_ins_vals)
}

