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
  evalue = double(), score = double(), bias = double()
)
for ( tbloutfile in opt$args ) {
  logmsg(sprintf("Reading %s", tbloutfile))
  tblout = tblout %>%
    union(
      read_fwf(
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
        ) %>% 
        select(accno, profile, evalue, score, bias)#,rest)
    )
}

logmsg("Done")
