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

# Args list for testing:
# opt = list(args = c('hmmsearch2classification.00.d/NrdAe.tblout','hmmsearch2classification.00.d/NrdAg.tblout','hmmsearch2classification.00.d/NrdAh.tblout','hmmsearch2classification.00.d/NrdAi.tblout','hmmsearch2classification.00.d/NrdAk.tblout','hmmsearch2classification.00.d/NrdAm.tblout','hmmsearch2classification.00.d/NrdAn.tblout','hmmsearch2classification.00.d/NrdAq.tblout','hmmsearch2classification.00.d/NrdA.tblout','hmmsearch2classification.00.d/NrdAz3.tblout','hmmsearch2classification.00.d/NrdAz4.tblout','hmmsearch2classification.00.d/NrdAz.tblout'), options=list(verbose=T))

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
  evalue = double(), score = double(), bias = double()
)
accessions = data.table(accno = character(), all = character())
for ( tbloutfile in opt$args ) {
  logmsg(sprintf("Reading %s", tbloutfile))
  t =  read_fwf(
    tbloutfile, fwf_cols(content = c(1, NA)), 
    col_types = cols(content = col_character()), 
    comment='#'
  ) %>% 
    separate(
      content, 
      c('accno', 't0', 'profile', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'), 
      '\\s+', 
      extra='merge',
      convert = T
    )
  tblout = union(tblout, t %>% select(accno, profile, evalue, score, bias))
  accessions = union(accessions, t %>% transmute(accno, all=sprintf("%s %s", accno, rest)))
}

# Make the accessions table a long map
accmap = data.table(accno=character(), accto=character())
while ( length((accessions %>% filter(!is.na(all)))$accno) > 0) {
  a = accessions %>% 
    separate(all, c(paste('c', 0:9), 'all'), sep='\x01', extra='merge', fill='right') %>% 
    gather(c, accto, 2:11) %>% 
    filter(!is.na(accto)) %>% 
    select(-c)
  accmap = union(accmap, a %>% select(-all))
  accessions = a %>% filter(!is.na(all)) %>% select(accno, all)
}

# For safety's sake: do a distinct on the accession map
accmap = accmap %>% distinct()

logmsg("Done")
