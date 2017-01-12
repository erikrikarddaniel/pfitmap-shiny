#!/usr/bin/env Rscript

# Converts all tsv files in the current directory to feather format.

library(readr)
library(feather)

for ( tsvfile in Sys.glob('*.tsv') ) {
  write(sprintf("Converting %s to feather format", tsvfile), stderr())
  write_feather(
    read_tsv(tsvfile),
    sub('\\.tsv', '.feather', tsvfile)
  )
}