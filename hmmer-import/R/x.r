domtblout <- read_tsv(
  'overlap_test.domtblout.tsv', 
  col_types = cols(.default = col_integer(), accno = col_character(), profile = col_character())
)

# Define a table that will be filled with lengths
lengths <- tibble(
  accno = character(), profile = character(), lentype = character(), len = integer()
)

# Function that joins all with n > 1 with the next row
do_nextjoin <- function(dt) {
  dt %>% 
    filter(n > 1) %>%
    left_join(
      domt %>% transmute(accno, profile, i = i - 1, next_row = TRUE, next_from = from, next_to = to), 
      by = c('accno', 'profile', 'i')
    ) %>%
    replace_na(list('next_row' = FALSE)) %>%
    return()
}

# FOR EACH FROM-TO PAIR:
for ( fs in list(
  c('hmm_from', 'hmm_to', 'hmmlen'),
  c('ali_from', 'ali_to', 'alilen'),
  c('env_from', 'env_to', 'envlen')
)) {
  logmsg(sprintf("Calculating %s", fs[3]))

  # 1. Set from and to to current pair, and calculate i (rownumber) and n (total domains) for each combination of accno and profile
  domt <- domtblout %>% transmute(accno, profile, from = .data[[fs[1]]], to = .data[[fs[2]]]) %>% 
    group_by(accno, profile) %>% mutate(n = n(), i = rank(from)) %>% ungroup()

  # 2. Move rows with n == 1 to nooverlaps
  nooverlaps <- domt %>% filter(n == 1) %>% select(accno, profile, from, to)

  while ( domt %>% filter(n > 1) %>% nrow() > 0 ) {
    logmsg(sprintf("Working on overlaps, nrow: %d", domt %>% nrow()))
    
    nextjoin <- do_nextjoin(domt)

    # 4. Move non-overlapping to nooverlaps
    nooverlaps <- nooverlaps %>%
      union(
        nextjoin %>%
          filter(next_from > to) %>%
          transmute(accno, profile, from = next_from, to = next_to)
      )

    # 5. Calculate a new domt by selecting the overlapping and moving the content from next_to to to
    domt <- nextjoin %>%
      filter(next_from <= to) %>%
      transmute(accno, profile, from, to = ifelse(next_to > to, next_to, to)) %>%
      group_by(accno, profile) %>% mutate(n = n(), i = rank(from)) %>% ungroup()
    
    nooverlaps <- nooverlaps %>%
      union(domt %>% filter(n == 1) %>% select(accno, profile, from, to))
  }
  
  lengths <- lengths %>%
    union(
      nooverlaps %>% mutate(len = to - from + 1) %>%
        group_by(accno, profile) %>% summarise(len = sum(len)) %>% ungroup() %>%
        mutate(lentype = fs[1])
    )
}

# Make the length table long so it can be joined with ...
lengths <- lengths %>% spread(lentype, len, fill = 0)