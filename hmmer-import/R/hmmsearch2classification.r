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
    c('--profilehierarchies', default=NA, help='tsv file with profile hiearchies')
  ),
  make_option(
    c('--singletable', default=NA, help='Write data in a single tsv format to this filename.')
  ),
  make_option(
    c('--taxflat', default=NA, help='Name of NCBI taxon table in "taxflat" format (see https://github.com/erikrikarddaniel/taxdata2taxflat).')
  ),
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
# opt = list(args = c('hmmsearch2classification.00.d/NrdAe.tblout','hmmsearch2classification.00.d/NrdAg.tblout','hmmsearch2classification.00.d/NrdAh.tblout','hmmsearch2classification.00.d/NrdAi.tblout','hmmsearch2classification.00.d/NrdAk.tblout','hmmsearch2classification.00.d/NrdAm.tblout','hmmsearch2classification.00.d/NrdAn.tblout','hmmsearch2classification.00.d/NrdAq.tblout','hmmsearch2classification.00.d/NrdA.tblout','hmmsearch2classification.00.d/NrdAz3.tblout','hmmsearch2classification.00.d/NrdAz4.tblout','hmmsearch2classification.00.d/NrdAz.tblout'), options=list(verbose=T, singletable='test.out.tsv', profilerhierarchies='/tmp/hmmtest.tsv', taxflat='taxflat.tsv'))

logmsg = function(msg, llevel='INFO') {
  if ( opt$options$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}
logmsg("Starting classification")

tblout = tibble(
  accno = character(), profile = character(),
  evalue = double(), score = double(), bias = double()
)
accessions = tibble(accno = character(), all = character())

# Read all the tblout files
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
accmap = data.table(accno=character(), accto=character(), desc=character(), taxon=character())
while ( length((accessions %>% filter(!is.na(all)))$accno) > 0) {
  a = accessions %>% 
    separate(all, c(paste('c', 0:9), 'all'), sep='\x01', extra='merge', fill='right') %>% 
    gather(c, accto, 2:11) %>% 
    filter(!is.na(accto)) %>% 
    select(-c)
  accmap = union(
    accmap, 
    a %>% select(-all) %>% 
      separate(accto, c('accto', 'desc'), sep=' ', extra='merge') %>%
      mutate(taxon=ifelse(
        grepl('\\[(.*)\\]', desc), 
        sub('.*\\[(.*)\\].*', '\\1', desc),
        'unknown')
      ) 
  )
  accessions = a %>% filter(!is.na(all)) %>% select(accno, all)
}

# For safety's sake: do a distinct on the accession map
accmap = accmap %>% distinct()

# Infer databases from the structure of accession numbers
accmap = accmap %>%
  mutate(db = ifelse(grepl('^.._', accto), 'refseq', NA)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[0-9A-Z]{4,4}_[0-9A-Z]$', accto)), 'pdb', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^P[0-9]+\\.[0-9]+$', accto)), 'uniprot', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[O,P,Q][0-9][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[A-NR-Z][0-9][A-Z][A-Z0-9][A-Z0-9][0-9]\\.[0-9]+$', accto)), 'uniprot', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[ADEKOJMNP][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'genbank', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[C][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'embl', db)) %>%
  mutate(db = ifelse((is.na(db) & grepl('^[BFGIL][A-Z][A-Z][0-9]+\\.[0-9]+$', accto)), 'dbj', db))

# Calculate best scoring profile for each accession
bestscoring = tblout %>% group_by(accno) %>% top_n(1, score) %>% ungroup()

# Join bestscoring with accmap to get a single table output
singletable = bestscoring %>% inner_join(accmap, by='accno') %>%
  transmute(db, accno=accto, profile, taxon, score, evalue)

# If we have a profile hierarchy file name, read it and join
if ( ! is.na(opt$options$profilehierarchies) ) {
  logmsg(sprintf("Adding profile hierarchies from %s", opt$options$profilehierarchies))
  singletable = singletable %>% left_join(
      read_tsv(opt$options$profilehierarchies, col_types=cols(.default=col_character())),
      by='profile'
    ) %>%
    select(-profile)
}

# If we have a taxflat NCBI taxonomy, read and join
if ( ! is.na(opt$options$taxflat) ) {
  logmsg(sprintf("Adding NCBI taxon ids from %s", opt$options$taxflat))
  singletable = singletable %>% 
    left_join(
      read_tsv(opt$options$taxflat, col_types=cols(.default=col_character(), ncbi_taxon_id=col_integer())) %>%
        select(taxon, ncbi_taxon_id),
      by='taxon'
    )
}

logmsg(sprintf("Writing single table %s", opt$options$singletable))
write_tsv(
  singletable %>% 
    arrange(accno),
  opt$options$singletable
)

logmsg("Done")
