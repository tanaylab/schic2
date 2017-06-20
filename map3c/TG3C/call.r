#! /net/mraid14/export/data/users/eladch/tools/CO6/R/3.2.1/bin/Rscript
#system("module update")

options(warn=1)

# get script name
all.args = commandArgs(F)
fn.arg = "--file="
script.name = sub(fn.arg, "", all.args[grep(fn.arg, all.args)])

args = commandArgs(T)
if (length(args) == 0) {
  cat(sprintf("usage: %s <working dir> <r file> <function> [parameters ...]\n", script.name))
  q(status=1) 
}

wd = args[1]
source.fn = args[2]
func = args[3]
params = args[4:length(args)]

# parse parameters
param.list = list()
for (i in seq_along(params))
{
  split = strsplit(params[i], "=")
  if (length(split[[1]]) != 2)
    stop(paste("value of", key, "is not defined"))
  key = split[[1]][1]
  value = split[[1]][2]
  
  if (is.character(value) && grepl(" ", value)) {
    value = strsplit(value, " ")[[1]]
  }

  if (all(is.element(value, c("T", "F"))))
    value = (value == "T")
  
  options(warn=-1) # to disable warning 'NAs introduced by coercion'
  if (all(!is.na(as.numeric(value))))
    value = as.numeric(value)
  options(warn=1)
  
  # cat(sprintf("key: %s, value: %s, type: %s\n", key, paste(value, collapse=","), typeof(value)))

  param.list[[key]] = value
}

# set working dir
setwd(wd)

# load source file
cat(sprintf("loading %s\n", source.fn))
suppressPackageStartupMessages(source(source.fn))
options(error=NULL)

tostr = function(x)
{
  if (typeof(x) == "character")
    x = paste("\"", x, "\"", sep="")
  if (length(x) > 1)
    x = paste("c(", paste(x, collapse=",", sep=""), ")", sep="")
  x
}

# print the function call, helpful for debug
fmt.list = lapply(param.list, tostr)
str = paste(names(fmt.list), fmt.list, sep="=", collapse=", ")
cat(sprintf("calling: %s(%s)\n", func, str))

# call the function
do.call(func, param.list)
