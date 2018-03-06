#!/usr/bin/env Rscript

# select_db_pclass_tdomain_spread_db_pool.r
#
# This script opens a pool of connections to an SQLite database, joins four tables, subsets on
# three columns (db, pclass and tdomain) and spreads the result.
#
# It expects the current directory to contain an sqlite file called: hmmsearch2classification.sqlite3.
#
# Author: daniel.lundin@dbb.su.se

suppressPackageStartupMessages(library(profvis))
suppressPackageStartupMessages(library(pool))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbplyr))
suppressPackageStartupMessages(library(tidyr))

dbpool <- dbPool(
  drv    = RSQLite::SQLite(),
  dbname = 'hmmsearch2classification.sqlite3'
)

# Setting up variables to connect to some tables
accessions           <- dbpool %>% tbl('accessions')
accessions_mem       <- accessions %>% collect()
hmm_profiles         <- dbpool %>% tbl('hmm_profiles')
hmm_profiles_mem     <- hmm_profiles %>% collect()
dupfree_proteins     <- dbpool %>% tbl('dupfree_proteins')
dupfree_proteins_mem <- dupfree_proteins %>% collect()
taxa                 <- dbpool %>% tbl('taxa')
taxa_mem             <- taxa %>% collect()

# Collect late
write("Collecting late", stderr())
system.time(
  {
  d <- dupfree_proteins %>%
    inner_join(accessions, by = 'accno') %>%
    inner_join(hmm_profiles, by = 'profile') %>%
    inner_join(taxa, by = 'taxon') %>%
    filter(db == 'refseq', pfamily == 'NrdJA', tdomain %in% c('Bacteria', 'Archaea', 'Eukaryota')) %>%
    group_by(pclass, tdomain) %>% summarise(n = n()) %>% ungroup() %>%
    collect() %>%
    spread(tdomain, n, fill = 0)
  }
)

# Collect after joins
write("Collecting after joins, before filter and group_by", stderr())
system.time(
  {
  d <- dupfree_proteins %>%
    inner_join(accessions, by = 'accno') %>%
    inner_join(hmm_profiles, by = 'profile') %>%
    inner_join(taxa, by = 'taxon') %>%
    collect() %>%
    filter(db == 'refseq', pfamily == 'NrdJA', tdomain %in% c('Bacteria', 'Archaea', 'Eukaryota')) %>%
    group_by(pclass, tdomain) %>% summarise(n = n()) %>% ungroup() %>%
    spread(tdomain, n, fill = 0)
  }
)

# Collect individual tables
write("Collecting individual tables", stderr())
system.time(
  {
  d <- dupfree_proteins %>% collect() %>%
    inner_join(accessions %>% collect(), by = 'accno') %>%
    inner_join(hmm_profiles %>% collect(), by = 'profile') %>%
    inner_join(taxa %>% collect(), by = 'taxon') %>%
    collect() %>%
    filter(db == 'refseq', pfamily == 'NrdJA', tdomain %in% c('Bacteria', 'Archaea', 'Eukaryota')) %>%
    group_by(pclass, tdomain) %>% summarise(n = n()) %>% ungroup() %>%
    spread(tdomain, n, fill = 0)
  }
)

# From in memory tables
write("From in memory tables", stderr())
system.time(
  {
  d <- dupfree_proteins_mem %>%
    inner_join(accessions_mem, by = 'accno') %>%
    inner_join(hmm_profiles_mem, by = 'profile') %>%
    inner_join(taxa_mem, by = 'taxon') %>%
    filter(db == 'refseq', pfamily == 'NrdJA', tdomain %in% c('Bacteria', 'Archaea', 'Eukaryota')) %>%
    group_by(pclass, tdomain) %>% summarise(n = n()) %>% ungroup() %>%
    collect() %>%
    spread(tdomain, n, fill = 0)
  }
)

poolClose(dbpool)
