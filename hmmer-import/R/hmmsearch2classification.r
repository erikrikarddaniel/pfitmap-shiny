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
    c('--profilehierarchies', default='', help='tsv file with profile hiearchies')
  ),
  make_option(
    c('--singletable', default='', help='Write data in a single tsv format to this filename.')
  ),
  make_option(
    c('--sqlitedb', default='', help='Write data in a SQLite database with this filename.')
  ),
  make_option(
    c('--taxflat', default='', help='Name of NCBI taxon table in "taxflat" format (see https://github.com/erikrikarddaniel/taxdata2taxflat).')
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

if ( length(grep('sqlitedb', names(opt$options), value = TRUE)) > 0 ) {
  suppressPackageStartupMessages(library(dbplyr))
}

# Args list for testing:
# opt = list(args = c('hmmsearch2classification.00.d/NrdAe.tblout','hmmsearch2classification.00.d/NrdAg.tblout','hmmsearch2classification.00.d/NrdAh.tblout','hmmsearch2classification.00.d/NrdAi.tblout','hmmsearch2classification.00.d/NrdAk.tblout','hmmsearch2classification.00.d/NrdAm.tblout','hmmsearch2classification.00.d/NrdAn.tblout','hmmsearch2classification.00.d/NrdAq.tblout','hmmsearch2classification.00.d/NrdA.tblout','hmmsearch2classification.00.d/NrdAz3.tblout','hmmsearch2classification.00.d/NrdAz4.tblout','hmmsearch2classification.00.d/NrdAz.tblout','hmmsearch2classification.00.d/NrdAe.domtblout','hmmsearch2classification.00.d/NrdAg.domtblout','hmmsearch2classification.00.d/NrdAh.domtblout','hmmsearch2classification.00.d/NrdAi.domtblout','hmmsearch2classification.00.d/NrdAk.domtblout','hmmsearch2classification.00.d/NrdAm.domtblout','hmmsearch2classification.00.d/NrdAn.domtblout','hmmsearch2classification.00.d/NrdAq.domtblout','hmmsearch2classification.00.d/NrdA.domtblout','hmmsearch2classification.00.d/NrdAz3.domtblout','hmmsearch2classification.00.d/NrdAz4.domtblout','hmmsearch2classification.00.d/NrdAz.domtblout'), options=list(verbose=T, singletable='test.out.tsv', profilehierarchies='hmmsearch2classification.00.phier.tsv', taxflat='taxflat.tsv', sqlitedb='testdb.sqlite3'))

logmsg = function(msg, llevel='INFO') {
  if ( opt$options$verbose ) {
    write(
      sprintf("%s: %s: %s", llevel, format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg),
      stderr()
    )
  }
}
logmsg("Starting classification")

# We will populate two tables, one with the full results, one with accessions
tblout = tibble(
  accno = character(), profile = character(),
  evalue = double(), score = double(), bias = double()
)
accessions = tibble(accno = character(), accto = character())

# Read all the tblout files
for ( tbloutfile in grep('\\.tblout', opt$args, value=TRUE) ) {
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
  accessions = union(accessions, t %>% transmute(accno, accto = sprintf("%s %s", accno, rest)))
}

# Split the accto field
accessions = accessions %>% separate_rows(accto, sep = '\x01') %>% 
  mutate(
    taxon = ifelse(
      grepl('\\[(.*)\\]', accto), 
      sub('.*\\[(.*)\\].*', '\\1', accto),
      'unknown'
    ),
    accto = sub(' .*', '', accto)
  )

domtblout = tibble(
  accno = character(), tlen = integer(), profile = character(), qlen = integer(), i = integer(), n = integer(), 
  dom_c_evalue = double(), dom_i_evalue = double(), dom_score = double(), dom_bias = double(),
  hmm_from = integer(), hmm_to = integer(), ali_from = integer(), ali_to = integer(), 
  env_from = integer(), env_to = integer()
)

# Read all the domtblout files
for ( domtbloutfile in grep('\\.domtblout', opt$args, value=TRUE) ) {
  logmsg(sprintf("Reading %s", domtbloutfile))
  t = read_fwf(
    domtbloutfile, fwf_cols(content = c(1, NA)), 
    col_types = cols(content = col_character()), 
    comment='#'
  ) %>% 
    separate(
      content, 
      c(
        'accno', 't0', 'tlen', 'profile', 't1', 'qlen',  'evalue', 'score', 'bias', 'i', 'n', 
        'dom_c_evalue', 'dom_i_evalue', 'dom_score', 'dom_bias', 
        'hmm_from', 'hmm_to', 'ali_from', 'ali_to', 'env_from', 'env_to', 'acc', 'rest'
      ),
      '\\s+', 
      extra='merge',
      convert = T
    )
  
  domtblout = union(
    domtblout,
    t %>% select(
      accno, tlen, profile, qlen, i, n, dom_c_evalue, dom_i_evalue, dom_score, dom_bias,
      hmm_from, hmm_to, ali_from, ali_to, env_from, env_to
    )
  )
}

# Calculate lengths by summing the ali_from and ali_to fields

# First check if there are overlaps...
domtblout.no_overlaps = domtblout %>%
  select(accno, profile, tlen, qlen, i, n, ali_from, ali_to)

calc_overlaps = function(dt) {
  o = dt %>% filter(n > 1) %>%
    select(accno, profile, ali_from, ali_to, i, n) %>%
    inner_join(
      domtblout %>% transmute(accno, profile, next_ali_from = ali_from, next_ali_to = ali_to, next_i = i, i = i - 1),
      by = c('accno', 'profile', 'i')
    ) %>% filter(next_ali_from <= ali_to) %>%
    mutate(new_ali_from = ali_from, new_ali_to = next_ali_to, new_n = n - 1)
  return(o)
}

overlaps = calc_overlaps(domtblout.no_overlaps)

while ( overlaps %>% nrow() > 0 ) {
  logmsg(sprintf("Handling overlaps, %d rows remaning", overlaps %>% nrow()))

  domtblout.no_overlaps = domtblout.no_overlaps %>% 
    # 1. Join in overlaps and set new ali_from and ali_to for matching rows
    left_join(
      overlaps %>% select(accno, profile, i, new_ali_from, new_ali_to), 
      by=c('accno', 'profile', 'i')
    ) %>%
    mutate(
      ali_from = ifelse(! is.na(new_ali_from), new_ali_from, ali_from),
      ali_to = ifelse(! is.na(new_ali_to), new_ali_to, ali_to)
    ) %>%
    # 2. Join in overlaps with i + 1 to get rid of the replaced row
    left_join(
      overlaps %>% transmute(accno, profile, i = next_i, new_n),
      by=c('accno', 'profile', 'i')
    ) %>%
    filter(is.na(new_n)) %>% select(-new_ali_from, -new_ali_to, -new_n) %>%
    # 3. Join in overlaps on only accno and profile to set new n for matching rows
    left_join(
      overlaps %>% select(accno, profile, new_n), by = c('accno', 'profile')
    ) %>%
    mutate(n = ifelse(!is.na(new_n), new_n, n)) %>% select(-new_n) %>%
    distinct(accno, profile, tlen, qlen, n, ali_from, ali_to) %>%
    # 4. Calculate new i
    group_by(accno, profile, tlen, qlen) %>% mutate(i = rank(ali_from)) %>% ungroup() %>%
    arrange(accno, profile, i)

  overlaps = calc_overlaps(domtblout.no_overlaps)
}

logmsg("Overlaps done")

# Now, we can calculate lengths
align_lengths = domtblout.no_overlaps %>%
  mutate(alilen = ali_to - ali_from + 1) %>%
  group_by(accno, profile, tlen, qlen) %>% summarise(alilen = sum(alilen)) %>% ungroup()

logmsg("Calculated lengths")

# Infer databases from the structure of accession numbers
accessions = accessions %>%
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
bestscoring = tblout %>% group_by(accno) %>% top_n(1, score) %>% ungroup() %>%
  select(accno, profile, score, evalue)

# Join in lengths
logmsg(sprintf("Joining in lengths from domtblout, nrows before: %d", bestscoring %>% nrow()))
bestscoring = bestscoring %>% inner_join(align_lengths, by = c('accno', 'profile'))

# If we were called with the singletable option, prepare data suitable for that
if ( length(grep('singletable', names(opt$options), value = TRUE)) > 0 ) {
  # If we have a profile hierarchy file name, read it and join
  if ( length(grep('profilehierarchies', names(opt$options), value = TRUE)) > 0 ) {
    logmsg(sprintf("Adding profile hierarchies from %s, nrows before: %d", opt$options$profilehierarchies, bestscoring %>% nrow()))
    bestscoring = bestscoring %>% 
      left_join(
        read_tsv(opt$options$profilehierarchies, col_types=cols(.default=col_character())),
        by='profile'
      )
  }

  # Join bestscoring with accessions and drop profile to get a single table output
  logmsg(sprintf("Joining in all accession numbers and dropping profile column, nrows before: %d", bestscoring %>% nrow()))
  singletable = bestscoring %>% inner_join(accessions, by='accno') %>%
    mutate(accno = accto) %>% select(-profile)

  # If we have a taxflat NCBI taxonomy, read and join
  if ( length(grep('taxflat', names(opt$options), value = TRUE)) > 0 ) {
    logmsg(sprintf("Adding NCBI taxon ids from %s, nrows before: %d", opt$options$taxflat, singletable %>% nrow()))
    singletable = singletable %>% 
      left_join(
        read_tsv(opt$options$taxflat, col_types=cols(.default=col_character(), ncbi_taxon_id=col_integer())) %>%
          select(taxon, ncbi_taxon_id),
        by='taxon'
      )
  }

  logmsg(sprintf("Writing single table %s, nrows: %d", opt$options$singletable, singletable %>% nrow()))
  write_tsv(
    singletable %>% 
      select(db, accno, taxon, score, evalue, psuperfamily:pgroup, ncbi_taxon_id, tlen:alilen) %>%
      arrange(accno),
    opt$options$singletable
  )
}

# If the user specified a filename for a SQLite database, write that here
if ( length(grep('sqlitedb', names(opt$options), value = TRUE)) > 0 ) {
  logmsg(sprintf("Creating/opening SQLite database %s", opt$options$sqlitedb))
  con = DBI::dbConnect(RSQLite::SQLite(), opt$options$sqlitedb, create = TRUE)
  
  con %>% copy_to(accessions, 'accessions', temporary = FALSE, overwrite = TRUE)
  
  con %>% copy_to(bestscoring, 'bestscoring', temporary = FALSE, overwrite = TRUE)

  if ( length(grep('profilehierarchies', names(opt$options), value = TRUE)) > 0 ) {
    logmsg(sprintf("Adding profile hierarchies from %s", opt$options$profilehierarchies))
    con %>% copy_to(
      read_tsv(opt$options$profilehierarchies, col_types=cols(.default=col_character())),
      'hmm_profiles', temporary = FALSE, overwrite = TRUE
    )
  }

  # If we have a taxflat NCBI taxonomy, read and join
  if ( length(grep('taxflat', names(opt$options), value = TRUE)) > 0 ) {
    logmsg(sprintf("Adding NCBI taxon ids from %s", opt$options$taxflat))
    con %>% copy_to(
      read_tsv(opt$options$taxflat, col_types=cols(.default=col_character(), ncbi_taxon_id=col_integer())) %>%
        transmute(
          ncbi_taxon_id,
          psuperkingdom = superkingdom, pkingdom = kingdom,
          pphylum       = phylum,       pclass   = class,
          porder        = order,        pfamily  = family,
          pgenus        = genus,        pspecies = species,
          taxon, rank
        ) %>%
        inner_join(accessions %>% distinct(taxon), by='taxon'),
      'taxa', temporary = FALSE, overwrite = TRUE
    )
  }
}

logmsg("Done")
