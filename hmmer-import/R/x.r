domtblout <- read_tsv(
  'overlap_test.domtblout.tsv', 
  col_types = cols(.default = col_integer(), accno = col_character(), profile = col_character())
)


# Join all with n > 1 with the next row
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

# 1. Set from and to to current pair, and calculate i (rownumber) and n (total domains) for each combination of accno and profile
domt <- domtblout %>% transmute(accno, profile, from = hmm_from, to = hmm_to) %>% 
  group_by(accno, profile) %>% mutate(n = n(), i = rank(from)) %>% ungroup()

# 2. Move rows with n == 1 to nooverlaps
nooverlaps <- domt %>% filter(n == 1) %>% select(accno, profile, from, to)

while ( domt %>% filter(n > 1) %>% nrow() > 0 ) {
  write(sprintf("domt > nrow(): %d", domt %>% nrow()), stderr())
  
  nextjoin <- do_nextjoin(domt)

  # 4. Move non-overlapping to nooverlaps
  nooverlaps <- nooverlaps %>%
    union(
      nextjoin %>%
        filter(next_from > to) %>%
        transmute(accno, profile, from = next_from, to = next_to)
    )
  write(sprintf("nooverlaps > nrow(): %d", nooverlaps %>% nrow()), stderr())

  # 5. Calculate a new domt by selecting the overlapping and moving the content from next_to to to
  domt <- nextjoin %>%
    filter(next_from <= to) %>%
    transmute(accno, profile, from, to = next_to) %>%
    group_by(accno, profile) %>% mutate(n = n(), i = rank(from)) %>% ungroup()
  
  nooverlaps <- nooverlaps %>%
    union(domt %>% filter(n == 1) %>% select(accno, profile, from, to))
}