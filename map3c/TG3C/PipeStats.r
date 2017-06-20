source(paste(Sys.getenv("PIPELINE_HOME"), "/map3c/TG3C/params.r", sep=""))

####
get.seg.stats = function(wd, seg.log.fn=NULL)
{
  if (is.null(seg.log.fn)) {
    seg.log.fn = paste(wd, "/logs/seg.log", sep="")
  }
  
  if (!file.exists(seg.log.fn)) {
    message("ERROR: cannot find log file ", seg.log.fn)
    return(NULL)
  }

  # parse counts from log line i this format: total:\t0: 16298010\t1: 1224283\t2: 2\t3: 
  lines = readLines(seg.log.fn, n=-1)
  total.line = grep(lines, pattern="^total:", value=T) 
  counts.str = strsplit(total.line, split="\t")[[1]][-1]

  counts = c()
  segs.in.read = c()
  for (x in counts.str) {
    x.s = strsplit(x, split=": ")[[1]]
    if (!is.na(x.s[2])) {
      segs.in.read = c(segs.in.read, x.s[1])
      counts = c(counts, as.numeric(x.s[2]))    
    }
  }
  names(counts) = segs.in.read
    
  counts
}

###### Set stats as track variables
update.tracks.stats=function(track.re.str, lanes.fn=NULL)
{
  tlist = gtrack.ls(track.re.str)

  if (!is.null(lanes.fn)) {
    lanes = read.delim(lanes.fn, header=F)
  }  

  for (tn in as.vector(tlist)) {
    cat(paste("Updating", tn, " stats..."))
    
    x = gextract(tn, gintervals=gintervals.2d.all())
    gtrack.var.set(tn, "n.adj", nrow(x)/2)
    gtrack.var.set(tn, "cis.frac", sum(as.character(x$chrom1) == as.character(x$chrom2))/nrow(x))

    if (!is.null(lanes.fn)) {
      value = as.character(lanes[which(lanes[,1] == tn), 2])
      if (length(value) > 1) {
        browser()
      }
      gtrack.var.set(tn, "lane", value)  
    }
    
    cat(" Done.\n")
  }
}

###### Get tracks stats
get.tracks.stats=function(track.re.str, vars=c("n.adj", "cis.frac", "lane", "early.cov.norm", "late.cov.norm"))
{
  tl = as.vector(gtrack.ls(track.re.str))
  res = data.frame(track=tl)
  for (v in vars) {
    res = cbind(res, sapply(tl, gtrack.var.get, var=v))
    colnames(res)[ncol(res)] = v
      
  }
  res
}

##################
detect.valid.4c.peaks.per.key=function(fid2key, min.fids.per.key=100, binsize=25e3, peak.bw=5, ofn="peaks.txt", frags.list="/net/mraid14/export/data/db/tgdb/hg19/seq/redb/GATC.frags", peak.factor=10, surround.factor=3, min.peaks.dist=5e5, mode="singleton", r.ofn="peaks.Rdata")
{
  chroms = gintervals.all()
  chroms$chrom = gsub(chroms$chrom, pattern="chr", replace="")
  
  cat ("reading tables ...\n")
  m = read.table(fid2key, header=T, stringsAsFactors=F)
  key.c = tapply(m$key, m$key, length)
  valid.keys = names(key.c[key.c >= min.fids.per.key])
  m = m[is.element(m$key, valid.keys), ]
  
  fids = read.delim(frags.list, header=F)
  colnames(fids) = c("fid", "chrom", "start", "end")
  fids$length = fids$end - fids$start
  
  cat ("grouping fids by KEY and chrom...\n")
  m = merge(m, fids)
  ms = split((m$start + m$end)/2, paste(m$key, m$chrom))

  peaks = NULL

  for (x in names(ms)) {
    key  = strsplit(x, split=" ")[[1]][1]
    chr = strsplit(x, split=" ")[[1]][2]
    v = ms[[x]]
    vt = table(ceiling(v / binsize) * binsize)

    cands = as.numeric(names(vt[vt >= peak.factor * median(vt)]))
    surround.mean = sapply(cands, function(x) { sum(vt[as.numeric(names(vt)) >= x - binsize * peak.bw & as.numeric(names(vt)) <= x + binsize * peak.bw & names(vt) != x])/(2 * peak.bw) })
    surround.mean[is.na(surround.mean)] = 0
      
    cands = cands[ surround.mean >= surround.factor * median(vt)]
    
    if (length(cands) > 0) {
      best.cands = c()
      cands = cands[order(vt[as.character(cands)], decreasing=T)]
      for (i in seq_along(cands)) {
        if (length(best.cands) == 0 || min(cands[i] - best.cands) >= min.peaks.dist) {
          best.cands = c(best.cands, cands[i])
        }
      }
      rep.fids = sapply(best.cands, select.rep.frag, fids=fids[fids[,2] == chr, ], binsize=binsize, exclude.fids=m[m$key == key & m$chrom == chr, 'fid'])
      
      peaks = rbind(peaks, data.frame(key=key, chrom=chr, coord=best.cands, frag=rep.fids, frag_start=fids[is.element(fids$fid, rep.fids), 'start'], frag_len=fids[is.element(fids$fid, rep.fids), 'length']))
    }
  }
  if (mode == "singleton") {
    key = table(peaks$key)
    s.key = names(key)[key == 1]
    peaks = peaks[is.element(peaks$key, s.key), ]
  }
  write.table(file=ofn, x=peaks, sep="\t", quote=F, row.names=F)
  if (!is.null(r.ofn)) {
    save(peaks, file=r.ofn)
  }
}


#####
get.frags.chrom.ranges=function(frags.fn, range.fn)
{
  if (!file.exists(range.fn)) {
    x = read.delim(frags.fn, header=F)
    from = tapply(x[,1], x[,2], min)
    to = tapply(x[,1], x[,2], max)
    r = data.frame(chrom=names(from), from=from, to=to)
    write.table(r, range.fn, quote=F, sep="\t")
  }
  read.table(range.fn, header=T, sep="\t")
}

select.rep.frag=function(coord, fids, binsize, exclude.fids)
{
  fids = fids[ !is.element(fids[,1], exclude.fids) & fids[,3] >= coord - binsize & fids[,4] <= coord, ]
  fids[ round(nrow(fids)/2), 1]
}
  
