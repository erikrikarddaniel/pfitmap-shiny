#!/usr/bin/env Rscript

library(data.table)
library(dtplyr)
library(dplyr)
library(readr)

cp = data.table(
  read_tsv(
    '~/data/pfitmap-eval/classified_proteins.prop_matching_ge_0.9.tsv',
    col_types = cols(
      .default = col_character(),
      profile_length = col_integer(),
      align_length = col_integer(),
      align_start = col_integer(),
      align_end = col_integer(),
      prop_matching = col_double(),
      ss_version = col_integer(),
      e_value = col_double(),
      score = col_double()
    )
  )
)

tcp = read_tsv(
  'classified_proteins.10000.tsv',
  col_types = cols(
    .default = col_character(),
    profile_length = col_integer(),
    align_length = col_integer(),
    align_start = col_integer(),
    align_end = col_integer(),
    prop_matching = col_double(),
    ss_version = col_integer(),
    e_value = col_double(),
    score = col_double()
  )
)

# Add Sulfolobus (run ONLY ONCE)
tcp = tcp %>%
  union(cp %>% filter(tspecies=='Sulfolobus acidocaldarius'))

write_tsv(tcp, 'classified_proteins.10000.tsv')

# This block would replace current test data with a new mostly random selection, DONT DO THAT.
### # Write a test file with ca. 10000 rows. The file will contain all proteins for a selection
### # of strains. The selection is random plus a few favourites.
### write_tsv(
###   cp %>% 
###     filter(db %in% c('ref', 'pdb', 'dbj')) %>%
###     inner_join(
###       cp %>% filter(db %in% c('ref', 'pdb', 'dbj')) %>% 
###         select(ncbi_taxon_id) %>% distinct() %>% sample_n(200) %>%
###         union(
###           cp %>% 
###             filter(db=='ref') %>% 
###             filter(
###               tspecies %in% c(
###                 'Escherichia coli', 'Pseudomonas aeruginosa', 
###                 'Thermotoga maritima', 'Enterobacteria phage T4',
###                 'Salmonella enterica', 'Lactobacillus leichmannii',
###                 'Clostridium botulinum',
###                 'Homo sapiens', 'Saccharomyces cerevisiae'
###               )
###             ) %>%
###             select(ncbi_taxon_id) %>% distinct()
###         ) %>%
###         union( tibble(ncbi_taxon_id = '382359') ) %>%
###         distinct(),
###       by = 'ncbi_taxon_id'
###     ),
###   'classified_proteins.10000.tsv'
### )

# Old, but keeping the code
### NrdGRE = cp %>% filter(psuperfamily=='NrdGRE')
### org = cp %>% select(tdomain:tstrain) %>% distinct()
### 
### NrdGRE.domclass = NrdGRE %>% group_by(tdomain, pclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain) %>% summarise(n_orgs=n()))
### write_tsv(NrdGRE.domclass, 'NrdGRE.domclass.tsv')
### 
### NrdGRE.domsubclass = NrdGRE %>% group_by(tdomain, pclass, psubclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain) %>% summarise(n_orgs=n()))
### write_tsv(NrdGRE.domsubclass, 'NrdGRE.domsubclass.tsv')
### 
### NrdGRE.orderclass = NrdGRE %>% group_by(tdomain,tkingdom,tphylum,tclass,torder, pclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain,tkingdom,tphylum,tclass,torder) %>% summarise(n_orgs=n()))
### write_tsv(NrdGRE.orderclass, 'NrdGRE.orderclass.tsv')
### 
### NrdGRE.ordersubclass = NrdGRE %>% group_by(tdomain,tkingdom,tphylum,tclass,torder, pclass, psubclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain,tkingdom,tphylum,tclass,torder) %>% summarise(n_orgs=n()))
### write_tsv(NrdGRE.ordersubclass, 'NrdGRE.ordersubclass.tsv')
### 
### NrdGREfer = cp %>% filter(psuperfamily %in% c('NrdGRE', 'Ferritin-like'))
### 
### NrdGREfer.domclass = NrdGREfer %>% group_by(tdomain, pclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain) %>% summarise(n_orgs=n()))
### write_tsv(NrdGREfer.domclass, 'NrdGREfer.domclass.tsv')
### 
### NrdGREfer.domsubclass = NrdGREfer %>% group_by(tdomain, pclass, psubclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain) %>% summarise(n_orgs=n()))
### write_tsv(NrdGREfer.domsubclass, 'NrdGREfer.domsubclass.tsv')
### 
### NrdGREfer.orderclass = NrdGREfer %>% group_by(tdomain,tkingdom,tphylum,tclass,torder, pclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain,tkingdom,tphylum,tclass,torder) %>% summarise(n_orgs=n()))
### write_tsv(NrdGREfer.orderclass, 'NrdGREfer.orderclass.tsv')
### 
### NrdGREfer.ordersubclass = NrdGREfer %>% group_by(tdomain,tkingdom,tphylum,tclass,torder, pclass, psubclass) %>% summarise(n_proteins = n()) %>% ungroup() %>% inner_join(org %>% group_by(tdomain,tkingdom,tphylum,tclass,torder) %>% summarise(n_orgs=n()))
### write_tsv(NrdGREfer.ordersubclass, 'NrdGREfer.ordersubclass.tsv')
