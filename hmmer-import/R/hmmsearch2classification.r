#!/usr/bin/env Rscript

# hmmsearch2classification
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dtplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# Get arguments
option_list = list(
  make_option(
    c("-v", "--verbose"), action="store_true", default=FALSE, 
    help="Print progress messages"
  )
)
opt = parse_args(
  OptionParser(option_list=option_list), 
  positional_arguments = TRUE
)

logmsg = function(msg, llevel='INFO') {
  if ( opt$options$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}
logmsg("Starting classification")

tblout = data.table(
  accno = character(), profile = character(),
  evalue = double(), score = double(), bias = double(),
  desc = character()
)
for ( tblout in opt$args ) {
  logmsg(sprintf("Reading %s", tblout))
  t = read_fwf('hmmsearch2classification.00.d/NrdAe.tblout', fwf_cols(content = c(1, NA)), comment='#') %>% separate(content, c('accno', 't0', 'profile', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'), '\\s+', extra='merge') %>% select(accno, profile, evalue, score, bias, rest)
}

logmsg("Done")
