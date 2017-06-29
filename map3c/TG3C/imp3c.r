source(paste(Sys.getenv("PIPELINE_HOME"), "/map3c/TG3C/params.r", sep=""))

#all paramters should be registered here
gtrack.create_mapab_track = function(params_fn)
{
  options(gparam.type = "string")
  Sys.setenv("PERL_BADLANG" = 0)
  params = init_params(params_fn)

  groot = get_param("TG3C.groot", params)
  
  gen_map_reads_pl = get_param("TG3C.gen_map_reads_pl", params)

  wd = get_param("TG3C.map_workdir", params)
  seq_fn = get_param("TG3C.full_seq_fn", params)
  
  rlen = as.numeric(get_param("TG3C.mapa_rlen", params))
  step = as.numeric(get_param("TG3C.mapa_step", params))
  binsize = as.numeric(get_param("TG3C.mapa_binsize", params))

  gdb.init(groot)
  
  #break the genome to reads
  message(sprintf("creating %d len reads from the genome\n", rlen))
  cmd = sprintf("perl %s %s %s %d %d -1", gen_map_reads_pl, seq_fn, wd, step, rlen)
  message(cmd)
  system(cmd)

  # map the reads to the genome
  message("invoking bowtie\n")
  bt2_bin = get_param("TG3C.bowtie2_bin", params)
  bt2_ndx = get_param("TG3C.bowtie2_ndx", params)
  nthreads = as.numeric(get_param("TG3C.bowtie2_threads", params))
	
  cmd = sprintf("%s -p %d --reorder --end-to-end -x %s -U %s/1.fastq -S %s/reads.sam 2>%s/bt.log", bt2_bin, nthreads, bt2_ndx, wd, wd, wd)
  message(cmd)
  system(cmd)

  # process mapped reads
  cmd = sprintf("cat %s/reads.sam | grep -v '^@' | awk 'OFS = \"\t\" { print $1, $3, $4, $10, \"+\"; } ' | sed 's/chr//g' | tr '_' '\t' | awk '$2 == $4 && $3 == $5' | cut -f 4-7 | sed 's/^/chr/g'  > %s/valid_reads", wd, wd, wd)
  message(cmd)
  system(cmd)
  
  # import to misha
  message("importing to misha ...")
  
  gdir.create("mapab", showWarnings=FALSE)

  ret <- gtrack.import_mappedseq(track=sprintf("mapab.length_%d", binsize),
                                 description=sprintf("mappability track rlen=%d step=%d binsize=%d", rlen, step, binsize),
                                 file=sprintf("%s/valid_reads", wd),
                                 pileup=binsize,
                                 binsize=step,
                                 cols.order=c(3,1,2,4)
                                 )
  message("Stats: ", paste(ret, collapse=" , "))

}

############
gtrack.create_redb_tracks = function(re_seq, params_fn, verbose=FALSE)
{
	options(gparam.type = "string")
	Sys.setenv("PERL_BADLANG" = 0)
	params = init_params(params_fn)

    groot = get_param("TG3C.groot", params)
    gdb.init(groot)
    
	gen_re_frags_pl = get_param("TG3C.gen_re_frags_pl", params)
	re_frags_to_fends_pl = get_param("TG3C.re_frags_to_fends_pl", params)
	mapab_track = get_param("TG3C.mapab_track", params)
	re_workdir = get_param("TG3C.re_workdir", params)
	chrom_key = get_param("TG3C.chrom_seq_key", params)

    gdir.create("redb")
    
	if(TG3C_miss_conf_err == TRUE) { 
		message("ERROR: missing configuration options\n")
		return(NA)
	}

#create the frags tables
	message("constructing table from genomic sequence\n")
	cmd = sprintf("perl %s %s %s %s", gen_re_frags_pl, re_seq, chrom_key, re_workdir)
    message(cmd)
	system(cmd)

#import the frags as sparse tracks
	message("importing fragments as sparse tracks")


	flen_track = sprintf("redb.%s_flen", re_seq)
	gc_track = sprintf("redb.%s_gc", re_seq)
	fragmap_track = sprintf("redb.%s_map", re_seq)
	flen_file = sprintf("%s/%s_flen", re_workdir, re_seq)
	gc_file = sprintf("%s/%s_gc", re_workdir, re_seq)
	fragmap_file = sprintf("%s/%s_map", re_workdir, re_seq)
	gtrack.import(flen_track, paste(re_seq, "fragment length"), flen_file, 0)
	gtrack.import(gc_track, paste(re_seq, "gc content"), gc_file, 0)

#annotate mapability
	message("generating mapability data for fraqments, using ", mapab_track)
	gtrack.create(fragmap_track, paste(re_seq, "fragment mapability"), 
		as.character(mapab_track), iterator=flen_track)

#dump mapability frags to text
	message("writing mapability text file to re workdir")
	gextract(fragmap_track, ALLGENOME, iterator=fragmap_track, file=fragmap_file)

#sperate frags to fends table
	message("writing fend tables from frag tables")
	cmd = sprintf("perl %s %s %s", re_frags_to_fends_pl, re_workdir, re_seq)
	system(cmd)

#import fend tables
	fe_flen_track = sprintf("redb.fe%s_flen", re_seq)
	fe_gc_track = sprintf("redb.fe%s_gc", re_seq)
	fe_fragmap_track = sprintf("redb.fe%s_map", re_seq)
	fe_flen_file = sprintf("%s/fe%s_flen", re_workdir, re_seq)
	fe_gc_file = sprintf("%s/fe%s_gc", re_workdir, re_seq)
	fe_fragmap_file = sprintf("%s/fe%s_map", re_workdir, re_seq)
	gtrack.import(fe_flen_track, paste(re_seq, "fragment end length"), fe_flen_file, 0)
	gtrack.import(fe_gc_track, paste(re_seq, "fragment end gc"), fe_gc_file, 0)
	gtrack.import(fe_fragmap_track,paste(re_seq, "fragment end mapability"), fe_fragmap_file, 0)
}


#########
#
gtrack.2d.create_from_3Cseq  = function(newtrack_name, params_fn, track_desc, groot=NULL, verbose=FALSE, combine_adj=FALSE, skip.if.exist=FALSE, overwrite.if.exist=TRUE, vars_dir="track_vars", interim_step=F)
{
  if (!is.null(groot)) {
    library(misha)
    gsetroot(groot)
  }
  message("1\n")
  options(gparam.type = "string")
  Sys.setenv("PERL_BADLANG" = 0)
  params = init_params(params_fn)
                                        #message(params)
  imp_3C_pipe_pl = get_param("TG3C.imp_3C_pipe_pl", params)
  merge_3C_pl = get_param("TG3C.merge_3C_pl", params)
  ndxs = get_param_list("TG3C.3C_indices", params)
  redb = get_param("TG3C.redb", params)
  re_seq = get_param("TG3C.RE_seq", params)
  fends_fn = sprintf("%s/%s.fends", redb, re_seq)

  work_dir = sprintf("%s/%s", get_param("TG3C.workdir", params), get_param("TG3C.3C_exp_nm", params));

  message("\n\nImport 3C from fastq\n")

  if(verbose) {
    message("work dir is ", work_dir)
  }

  if(TG3C_miss_conf_err == TRUE) { 
    message("ERROR: missing configuration options\n")
    return(NA)
  }

  if(!file.exists(fends_fn)) {
    message("ERROR: fends table is missing from ", fends_fn)
  }

  adj_list = c()
  for (foc_ndx in unlist(ndxs)) {
    if(verbose) { message("working on index ", foc_ndx, "\n") }
    
#    foc_ndx_track_name = sprintf("%s_%s", newtrack_name, foc_ndx)
    foc_ndx_track_name = sprintf("%s.%s_%s", newtrack_name, get_param("TG3C.3C_exp_nm", params), foc_ndx) # Y030216
    if (gtrack.exists(foc_ndx_track_name)) {
      if (skip.if.exist) {
        if (verbose) { message("skipping processing, track ", foc_ndx_track_name, " already exist.\n") }
        next
      }
      else if (overwrite.if.exist) {
        message("track ", foc_ndx_track_name, " already exist, overwriting it...\n")
        gtrack.rm(foc_ndx_track_name, T)
      }
      
    }
    
    cmd = sprintf("perl %s @%s -TG3C.foc_ndx %s\n", imp_3C_pipe_pl, params_fn, foc_ndx)
    if(verbose) { message(cmd,"\n") }
    system(cmd)
    
    status_fn = sprintf("%s.%s/done_ok", work_dir, foc_ndx)
    if(!interim_step && !file.exists(status_fn)) {
      message("ERROR: failed in index ", foc_ndx, "\nsee more logs in work dir ", work_dir)
      return(NA)
    }

                                        # import a track per foc_ndx adj file
    if (!interim_step && !combine_adj) { 
      contacts_fn = sprintf("%s.%s/adj", work_dir, foc_ndx)
      
      message("will import ", foc_ndx_track_name, " from ", contacts_fn, " fends at ", fends_fn)
      gtrack.2d.import_contacts(foc_ndx_track_name, description=track_desc, contacts=contacts_fn, fends=fends_fn)
      
      gtrack.import.vars.from.files(foc_ndx_track_name, sprintf("%s.%s/%s/", work_dir, foc_ndx, vars_dir))

      
    }
    else {
      adj_list = c(adj_list, sprintf("%s.%s/adj", work_dir, foc_ndx))
    }
                                        #compute summary stats from all ndxs
  }
  if (!interim_step && combine_adj) {
                                        #now combine all adj files
    message("will import ", newtrack_name, " from ", paste(adj_list, collape=", "), " fends at ", fends_fn)
    
    gtrack.2d.import_contacts(newtrack_name,
                              description=track_desc, 
                              contacts=adj_list,
                              fends=fends_fn)

    #gtrack.import.vars.from.files(newtrack_name, list.files(path=work_dir, pattern=vars.fn.regexp, full.names=T))
    
    return(newtrack_name)
  }
}


gtrack.import.vars.from.files = function(tn, idir)
{
  Sys.sleep(10)
  gdb.reload()
  r.ifns = list.files(path=idir, pattern=".Rdata", full.names=T)
  for (f in r.ifns) {
    vs = load(f)
    for (v in vs) {
      message("loading var ", v, " to track ", tn)
      gtrack.var.set(tn, v, get(v))
    }
  }

  tags.ifns = list.files(path=idir, pattern=".Rtable", full.names=T)
  for (f in tags.ifns) {
    tab = read.table(f, header=T)
    for (i in 1:nrow(tab)) {
      message("loading var ", tab[i, 'tag'], " to track ", tn)
      gtrack.var.set(tn, tab[i, 'tag'], tab[i, 'value'])
    }
  }
}

