TG3C_miss_conf_err = FALSE

init_params = function(fn)
{
  t = read.table(fn, sep="=", fill=T)

  params = sapply(as.character(t[,2]), FUN=function(v) { system2("echo", v, stdout=T) })
  names(params) = t[,1]
  assign("TG3C_miss_conf_err", FALSE, envir=.GlobalEnv)

  lines = readLines(fn, n=-1)
  inc.l = grep(lines, patt="^#include ", value=T)
  for (i in seq_along(inc.l)) {
    inc.fn = strsplit(inc.l[i], split=" +", perl=T)[[1]][2]
    if (file.exists(inc.fn)) {
      params = c(params, init_params(inc.fn))
    }
    else {
      message("cannot find file ", inc.fn, " to include")
    }
  }

  return(params)
}

get_param = function(nm, params) {

	if(nm %in% names(params)) {
		return(params[nm])
	} else {
		assign("TG3C_miss_conf_err", TRUE, envir=.GlobalEnv)
		message("missing params ", nm)
		return(NA)
	}
}

get_param_list= function(nm, params) {

	if(nm %in% names(params)) {
		return(strsplit(as.character(params[nm]),","))
	} else {
		assign("TG3C_miss_conf_err", TRUE, envir=.GlobalEnv)
		message("missing params ", nm)
		return(NA)
      }
}
