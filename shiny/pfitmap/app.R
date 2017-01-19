#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(dplyr)
library(data.table)
library(dtplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(stringr)
library(DT)
library(chorddiag)
library(feather)

# Some constants
PROTEIN_HIERARCHY = c( 'psuperfamily', 'pfamily', 'pclass', 'psubclass', 'pgroup' )
TAXON_HIERARCHY = c( 'tdomain', 'tkingdom', 'tphylum', 'tclass', 'torder', 'tfamily', 'tgenus', 'tspecies', 'tstrain' )

INDPROTEINS = 'indproteins'
COMBPROTEINS = 'combproteins'

DIV_PALETTE = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
DIV_PALETTE_96X = c(
  DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE, DIV_PALETTE
)
DIV_PALETTE_768X = c(
  DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X, DIV_PALETTE_96X
)
DIV_PALETTE_6144X = c(
  DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X, DIV_PALETTE_768X
)

LIGHT_PALETTE = c('#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99')
LIGHT_PALETTE_96X = c(
  LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE,
  LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE, LIGHT_PALETTE
)
LIGHT_PALETTE_768X = c(
  LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X, LIGHT_PALETTE_96X
)
LIGHT_PALETTE_6144X = c(
  LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X, LIGHT_PALETTE_768X
)
LIGHT_PALETTE_98304X = c(
  LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X,
  LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X, LIGHT_PALETTE_6144X
)

DARK_PALETTE = c('#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a','#b15928')
DARK_PALETTE_96X = c(
  DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE,
  DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE, DARK_PALETTE
)
DARK_PALETTE_768X = c(
  DARK_PALETTE_96X, DARK_PALETTE_96X, DARK_PALETTE_96X, DARK_PALETTE_96X, DARK_PALETTE_96X, DARK_PALETTE_96X, DARK_PALETTE_96X, DARK_PALETTE_96X
)

# Reading data and transforming
if ( grepl('\\.tsv$', Sys.getenv('PFITMAP_DATA')) ) {
  write(sprintf("LOG: %s: Reading tsv data from %s", Sys.time(), Sys.getenv('PFITMAP_DATA')), stderr())
  classified_proteins = data.table(
    read_tsv(
      Sys.getenv('PFITMAP_DATA'),
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
} else {
  if ( grepl('\\.feather$', Sys.getenv('PFITMAP_DATA')) ) {
    write(sprintf("LOG: %s: Reading feather data from %s", Sys.time(), Sys.getenv('PFITMAP_DATA')), stderr())
    classified_proteins = data.table(
      read_feather(Sys.getenv('PFITMAP_DATA'))
    )
  }
}

write(sprintf("LOG: %s: Filling in protein hierarchy", Sys.time()), stderr())
classified_proteins = data.table(
  classified_proteins %>%
    mutate(
      pfamily = ifelse(is.na(pfamily), sprintf("%s, no family", psuperfamily), pfamily),
      pclass = ifelse(is.na(pclass), sprintf("%s, no class", pfamily), pclass),
      psubclass = ifelse(is.na(psubclass), sprintf("%s, no subclass", pclass), psubclass),
      pgroup = ifelse(is.na(pgroup), sprintf("%s, no group", psubclass), pgroup)
    )
)

# We have problematic organisms, where multiple sequences of the same kind are
# assigned to the same taxon, a species or a genus. Trying to get rid of the
# species cases by filtering non-strain taxa from species with at least one
# strain. At the same time, delete all taxa with tgenus == tstrain and no
# tspecies.

write(sprintf("LOG: %s: Finding correct taxa", Sys.time()), stderr())

# Step 1. Get all unique taxa
taxa = data.table(
  classified_proteins %>% 
    select(db, ncbi_taxon_id, tdomain, tkingdom, tphylum, tclass, torder, tfamily, tgenus, tspecies, tstrain) %>% 
    distinct() %>% 
    filter( ! ( tgenus == tstrain & is.na(tspecies) ) ) %>%
    mutate(tspecies = ifelse(is.na(tspecies) & ! is.na(tgenus), sprintf("%s sp.", tgenus), tspecies))
)

# Step 2. Left join with a list of species that have strains, and then filter.
taxa = data.table(
  taxa %>%
    left_join(
      taxa %>% 
        filter(tspecies != tstrain) %>% 
        select(db,tdomain:tspecies) %>% distinct() %>% 
        mutate(strains=T),
      by = c("db", "tdomain", "tkingdom", "tphylum", "tclass", "torder", "tfamily", "tgenus", "tspecies")
    ) %>% 
    replace_na(list('strains'=F)) %>%
    filter( ! ( strains & tspecies == tstrain ) )
)

# Fill in empty levels of the taxon hierarchy (can't be done before the steps
# involving taxa above).
write(sprintf("LOG: %s: Filling in empty taxa in classified_proteins table", Sys.time()), stderr())
classified_proteins = data.table(
  classified_proteins %>%
    mutate(
      tkingdom = ifelse(is.na(tkingdom), sprintf("%s, no kingdom", tdomain), tkingdom),
      tphylum = ifelse(is.na(tphylum), sprintf("%s, no phylum", sub(', no kingdom', '', tkingdom)), tphylum),
      tclass = ifelse(is.na(tclass), sprintf("%s, no class", sub(', no phylum', '', tphylum)), tclass),
      torder = ifelse(is.na(torder), sprintf("%s, no order", sub(', no class', '', tclass)), torder),
      tfamily = ifelse(is.na(tfamily), sprintf("%s, no family", sub(', no order', '', torder)), tfamily),
      tgenus = ifelse(is.na(tgenus), sprintf("%s, no genus", sub(', no family', '', tfamily)), tgenus),
      tspecies = ifelse(is.na(tspecies), sprintf("%s, no species", sub(', no genus', '', tgenus)), tspecies)
    )
)

# Do the same for taxa
write(sprintf("LOG: %s: Filling in empty taxa in taxa table", Sys.time()), stderr())
taxa = data.table(
  taxa %>%
    mutate(
      tkingdom = ifelse(is.na(tkingdom), sprintf("%s, no kingdom", tdomain), tkingdom),
      tphylum = ifelse(is.na(tphylum), sprintf("%s, no phylum", sub(', no kingdom', '', tkingdom)), tphylum),
      tclass = ifelse(is.na(tclass), sprintf("%s, no class", sub(', no phylum', '', tphylum)), tclass),
      torder = ifelse(is.na(torder), sprintf("%s, no order", sub(', no class', '', tclass)), torder),
      tfamily = ifelse(is.na(tfamily), sprintf("%s, no family", sub(', no order', '', torder)), tfamily),
      tgenus = ifelse(is.na(tgenus), sprintf("%s, no genus", sub(', no family', '', tfamily)), tgenus),
      tspecies = ifelse(is.na(tspecies), sprintf("%s, no species", sub(', no genus', '', tgenus)), tspecies)
    )
)

# We will need a vector of protein superfamilies
psuperfamilies = (classified_proteins %>% select(psuperfamily) %>% distinct() %>% arrange(psuperfamily))$psuperfamily
tdomains = (taxa %>% select(tdomain) %>% distinct() %>% arrange(tdomain))$tdomain

# We also need a vector of databases
dbs = (classified_proteins %>% select(db) %>% distinct() %>% arrange(db))$db

write(sprintf("LOG: %s: Data init done", Sys.time()), stderr())

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel('pfitmap/RNRdb'),
  
  sidebarLayout(
    sidebarPanel( 
      radioButtons(
        'protstattype', 'Type of protein statistic',
        list(
          'Individual proteins' = INDPROTEINS,
          'Combinations of proteins' = COMBPROTEINS
        ),
        selected = INDPROTEINS
      ),
      selectInput(
        'db', 'Database',
        dbs, selected = 'ref'
      ),
      selectInput(
        'taxonrank', 'Taxon rank', 
        list(
          'Domain'  = 'tdomain',  'Phylum'  = 'tphylum', 'Class'   = 'tclass',
          'Order'   = 'torder',   'Family'  = 'tfamily', 'Genus'   = 'tgenus',
          'Species' = 'tspecies', 'Strain'  = 'tstrain'
        )
      ),
      selectInput(
        'proteinrank', 'Protein rank',
        list(
          'Superfamily' = 'psuperfamily', 'Family'   = 'pfamily',
          'Class'       = 'pclass',       'Subclass' = 'psubclass',
          'Group'       = 'pgroup'
        ),
        selected = c('pclass')
      ),
      wellPanel(
        selectInput(
          'psuperfamilies', 'Protein superfamilies',
          c('', psuperfamilies), selected = c('NrdGRE'), 
          multiple=T
        ),
        uiOutput('pfamilies'),
        uiOutput('pclasses')
      ),
      wellPanel(
        selectInput(
          'tdomains', 'Taxonomic domains',
          c('', tdomains), multiple = T
        ),
        uiOutput('tphyla'),
        uiOutput('tclasses'),
        uiOutput('torders'),
        uiOutput('tfamilies'),
        uiOutput('tgenera'),
        uiOutput('tspecies')
      )
    ),
    mainPanel(
      textOutput('ssversion'),
      textOutput('debug'),
      tabsetPanel(type= 'tabs', 
        tabPanel('table', 
          fluidRow(
            column(4, checkboxInput('taxonomysort', 'Taxonomic sort', value=T)),
            column(4, 
              selectInput(
                'trank4colour', 'Colour by taxon', 
                list('Domain' = 'tdomain'), selected='tdomain'
              )
            )
          ),
          dataTableOutput('mainmatrix')
        ),
        tabPanel('chord graph',
          chorddiagOutput('chordgraph', height=800)
        ),
        tabPanel('distributions',
          selectInput(
            'sinastat', 'Statistic',
            list(
              'HMM score' = 'score', 'Sequence length' = 'seqlen', 'Alignment length' = 'align_length'
            )
          ),
          plotOutput('distgraph', height=600)
        ),
        tabPanel('sequences',
          downloadLink('fastaseq', 'Download sequences in fasta format')
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # Reactive methods outside assignments
  
  # Returns a table after applying all filters the user have called for
  filtered_table = reactive({
    t = classified_proteins %>% 
      filter(db == input$db) %>%
      inner_join(taxa %>% select(db, ncbi_taxon_id), by=c('db', 'ncbi_taxon_id'))
    
    # Filters for protein hierarchy
    if ( length(input$psuperfamilies) > 0 ) { t = t %>% filter(psuperfamily %in% input$psuperfamilies) }
    if ( length(input$pfamilies) > 0 ) { t = t %>% filter(pfamily %in% input$pfamilies) }
    if ( length(input$pclasses) > 0 ) { t = t %>% filter(pclass %in% input$pclasses) }

    # Filters for taxon hierarchy
    if ( length(input$tdomains) > 0 ) { t = t %>% filter(tdomain %in% input$tdomains) }
    if ( length(input$tphyla) > 0 ) { t = t %>% filter(tphylum %in% input$tphyla) }
    if ( length(input$tclasses) > 0 ) { t = t %>% filter(tclass %in% input$tclasses) }
    if ( length(input$torders) > 0 ) { t = t %>% filter(torder %in% input$torders) }
    if ( length(input$tfamilies) > 0 ) { t = t %>% filter(tfamily %in% input$tfamilies) }
    if ( length(input$tgenera) > 0 ) { t = t %>% filter(tgenus %in% input$tgenera) }
    if ( length(input$tspecies) > 0 ) { t = t %>% filter(tspecies %in% input$tspecies) }

    # Construct a field for taxonomical sort and on for tooltip for taxonomy.
    # This is done in two steps: first a string is constructed, then the the
    # string is used in a mutate_ statement.
    n = which(TAXON_HIERARCHY==input$taxonrank)
    ts_string = paste(
      'sprintf(strrep("%-50s",', which(TAXON_HIERARCHY == input$taxonrank),  '), ',
      paste(TAXON_HIERARCHY[1:which(TAXON_HIERARCHY==input$taxonrank)], collapse=", "), ')'
    )
    ttt_string = paste(
      'paste(', paste(TAXON_HIERARCHY[1:which(TAXON_HIERARCHY==input$taxonrank)], collapse=", "), ', sep="; ")'
    )
    t = t %>%
      mutate_(
        'tsort' = ts_string, 
        'tcolour' = input$trank4colour,
        'taxon_tooltip' = ttt_string
      )
    
    write(sprintf("DEBUG: %s: rows in table: %d", Sys.time(), length(t[,1])), stderr())
    
    t
  })

  # Returns a filtered and summarised table after applying the group by
  # criteria called for by the user.
  indproteins_sums_table = reactive({
    d = filtered_table() %>%
      group_by_('tsort', 'tcolour', 'taxon_tooltip', input$taxonrank, input$proteinrank) %>%
      summarise(n=n()) %>%
      inner_join(
        taxa %>%
          inner_join(
            classified_proteins %>%
              filter(db == input$db) %>%
              select(db, ncbi_taxon_id) %>% distinct(),
            by = c('db', 'ncbi_taxon_id')
          ) %>%
          group_by_(input$taxonrank) %>%
          summarise(n_genomes = n()),
        by=c(input$taxonrank)
      ) %>%
      mutate(fraction = n/n_genomes) %>%
      gather(s, v, n, fraction) %>%
      mutate_(
        'proteinrank' = paste("ifelse(s == 'n',", input$proteinrank, ', sprintf("%s%s",', input$proteinrank, ', s))'),
        'n' = 'v'
      ) %>%
      select_(paste('-', input$proteinrank)) %>% select(-s, -v)
    
    write(sprintf("DEBUG: %s: rows in grouped table: %d", Sys.time(), length(d[,1])), stderr())
    
    d
  })

  # Calculates all present combinations of proteins at the proteinrank selected,
  # groups by the combinations, combines with n_genomes from taxa and returns.
  combproteins_sums_table = reactive({
    ###write(sprintf("combproteins_sums_table, protstattype: %s", input$protstattype), stderr())
    d = filtered_table() %>%
      group_by_('tsort', 'tcolour', 'taxon_tooltip', 'ncbi_taxon_id', input$taxonrank, input$proteinrank) %>%
      summarise(n=n()) %>%
      ungroup() %>%
      spread_(input$proteinrank, input$proteinrank, fill='') %>%
      select(-n)
    d = d %>%
      unite(comb, 6:length(colnames(d)), sep=':') %>%
      mutate(comb=sub(':::*', ':', sub(':*$', '', sub('^:*', '', comb)))) %>%
      group_by_('tsort', 'tcolour', 'taxon_tooltip', input$taxonrank, 'comb') %>%
      summarise(n=n()) %>%
      inner_join(
        taxa %>%
          inner_join(
            classified_proteins %>%
              filter(db == input$db) %>%
              select(db, ncbi_taxon_id) %>% distinct(),
            by = c('db', 'ncbi_taxon_id')
          ) %>%
          group_by_(input$taxonrank) %>%
          summarise(n_genomes = n()),
        by=c(input$taxonrank)
      ) %>%
      mutate(fraction = n/n_genomes) %>%
      gather(s, v, n, fraction) %>%
      mutate(
        comb = ifelse(s == 'n', comb, sprintf("%s%s", comb, s)),
        n = v
      ) %>%
      select(-s, -v)

    d
  })

  output$pfamilies = renderUI({
    pf = classified_proteins
    if ( length(input$psuperfamilies) > 0 ) {
      pf = pf %>% filter(psuperfamily %in% input$psuperfamilies)
    }
    pf = (pf %>% select(pfamily) %>% distinct())$pfamily
    
    selectInput(
      'pfamilies', 'Families',
      pf, multiple = T
    )
  })
  
  output$pclasses = renderUI({
    pc = classified_proteins
    if ( length(input$pfamilies) > 0) {
      pc = pc %>% filter(pfamily %in% input$pfamilies)
    }
    pc = (pc %>% select(pclass) %>% distinct() %>% arrange(pclass))$pclass
    
    selectInput(
      'pclasses', 'Classes',
      pc, multiple = T
    )
  })

  output$tphyla = renderUI({
    ###write(sprintf("input$tdomains, len: %d: %s", length(input$tdomains), input$tdomains), stderr())
    if ( length(input$tdomains) > 0 ) {
      selectInput(
        'tphyla', 'Phyla',
        (taxa %>% filter(tdomain %in% input$tdomains) %>% select(tphylum) %>% distinct() %>% arrange(tphylum))$tphylum,
        multiple = T
      )
    } else {
      selectInput(
        'tphyla', 'Phyla',
        (taxa %>% select(tphylum) %>% distinct() %>% arrange(tphylum))$tphylum,
        multiple = T
      )
    }
  })
  
  output$tclasses = renderUI({
    ###write(sprintf("tphyla: %d: %s", length(input$tphyla), input$tphyla), stderr())
    if ( length(input$tphyla) > 0 ) {
      selectInput(
        'tclasses', 'Classes',
        (taxa %>% filter(tphylum %in% input$tphyla) %>% select(tclass) %>% distinct() %>% arrange(tclass))$tclass,
        multiple = T
      )
    } else {
      selectInput(
        'tclasses', 'Classes',
        (taxa %>% select(tclass) %>% distinct() %>% arrange(tclass))$tclass,
        multiple = T
      )
    }
  })
  
  output$torders = renderUI({
    ###write(sprintf("tclasses: %d: %s", length(input$tclasses), input$tclasses), stderr())
    if ( length(input$tclasses) > 0 ) {
      selectInput(
        'torders', 'Orders',
        (taxa %>% filter(tclass %in% input$tclasses) %>% select(torder) %>% distinct() %>% arrange(torder))$torder,
        multiple = T
      )
    } else {
      selectInput(
        'torders', 'Orders',
        (taxa %>% select(torder) %>% distinct() %>% arrange(torder))$torder,
        multiple = T
      )
    }
  })
  
  output$tfamilies = renderUI({
    ###write(sprintf("torders: %d: %s", length(input$torders), input$torders), stderr())
    if ( length(input$torders) > 0 ) {
      selectInput(
        'tfamilies', 'Families',
        (taxa %>% filter(torder %in% input$torders) %>% select(tfamily) %>% distinct() %>% arrange(tfamily))$tfamily,
        multiple = T
      )
    } else {
      selectInput(
        'tfamilies', 'Families',
        (taxa %>% select(tfamily) %>% distinct() %>% arrange(tfamily))$tfamily,
        multiple = T
      )
    }
  })
  
  output$tgenera = renderUI({
    ###write(sprintf("tfamilies: %d: %s", length(input$tfamilies), input$tfamilies), stderr())
    if ( length(input$tfamilies) > 0 ) {
      selectInput(
        'tgenera', 'Genera',
        (taxa %>% filter(tfamily %in% input$tfamilies) %>% select(tgenus) %>% distinct() %>% arrange(tgenus))$tgenus,
        multiple = T
      )
    } else {
      selectInput(
        'tgenera', 'Genera',
        (taxa %>% select(tgenus) %>% distinct() %>% arrange(tgenus))$tgenus,
        multiple = T
      )
    }
  })
  
  output$tspecies = renderUI({
    if ( length(input$tgenera) > 0 ) {
      selectInput(
        'tspecies', 'Species',
        (taxa %>% filter(tgenus %in% input$tgenera) %>% select(tspecies) %>% distinct() %>% arrange(tspecies))$tspecies,
        multiple = T
      )
    } else {
      selectInput(
        'tspecies', 'Species',
        (taxa %>% select(tspecies) %>% distinct() %>% arrange(tspecies))$tspecies,
        multiple = T
      )
    }
  })

  output$trank4colour = renderUI({
    ranks = list()
    for ( r in TAXON_HIERARCHY[1:which(TAXON_HIERARCHY==input$taxonrank)] ) { 
      ranks[[Hmisc::capitalize(sub('^t', '', r))]] = r 
    }
    updateSelectInput(
      'trank4colour', 'Colour by taxon',
      ranks, selected = 'tdomain'
    )
  })
  
  output$mainmatrix = renderDataTable(
    {
      t = switch(
        input$protstattype,
        indproteins  = indproteins_sums_table() %>% spread(proteinrank, n, fill=0),
        combproteins = combproteins_sums_table() %>% spread(comb, n, fill=0)
      )
      write(sprintf("DEBUG: %s: \tcolnames: %s", Sys.time(), paste(colnames(t), collapse=", ")), stderr())
      #write.csv(t, stderr(), row.names=F)
      
      if ( input$taxonomysort ) {
        t = t %>% arrange(tsort)
      }
      
      # This is to get the right column names, a bit involved perhaps...
      t = t %>% mutate_('Taxon'=input$taxonrank, `N. genomes`='n_genomes') %>%
        mutate(Taxon = sprintf("<span title='%s'>%s</span>", taxon_tooltip, Taxon))
      c = colnames(t)
      write(sprintf("DEBUG: %s: \tcolnames: %s", Sys.time(), paste(colnames(t), collapse=", ")), stderr())
      
      # Colours for heatmap
      brks = quantile(c(0,1), probs = seq(.0, 1, .05), na.rm = TRUE)
      clrs = round(seq(255, 40, length.out = length(brks) + 1), 0) %>% { paste0("rgb(255,", ., ",", ., ")") }
      
      # Hide some columns
      invisible = c(0, grep('fraction', c) - 3)
      write(sprintf("--> invisible: %s", paste(invisible, collapse=", ")), stderr())

      write(sprintf("DEBUG: %s: \tcolnames: %s", Sys.time(), paste(colnames(t), collapse=", ")), stderr())
      t = t %>%
        select(tcolour, c(length(c)-1,length(c),7:length(c)-2))
      write(sprintf("DEBUG: %s: \tcolnames: %s", Sys.time(), paste(colnames(t), collapse=", ")), stderr())
      #write.csv(t, stderr(), row.names=F)
      dt = datatable(
        t, 
        rownames=F, 
        escape = c(T, F),
        options=list(
          lengthMenu = c(50, 100, 250, 500),
          columnDefs = list(list(targets = invisible, visible = FALSE))
        )
      ) %>%
        formatStyle(
          'Taxon', 'tcolour',
          backgroundColor = styleEqual(
            unique(t$tcolour), LIGHT_PALETTE_98304X[1:length(unique(t$tcolour))]
          )
        )
      
      for (ff in c[grep('fraction', c)]) {
        dt = dt %>% formatStyle(
          sub('fraction', '', ff), ff,
          backgroundColor = styleInterval(brks, clrs)
        )
      }
      
      dt
    }
  )
  
  output$chordgraph = renderChorddiag({
    t = switch(
      input$protstattype,
      indproteins  = indproteins_sums_table() %>% spread(proteinrank, n, fill=0),
      combproteins = combproteins_sums_table() %>% spread(comb, n, fill=0)
    ) %>% 
      select(-tsort, -tcolour, -taxon_tooltip, -n_genomes)
    
    # Delete all columns ending in 'fraction'
    for ( c in colnames(t)[grep('fraction', colnames(t))] ) {
      t = t %>% select_(sprintf("-`%s`", c))
    }
    
    ###write(sprintf("colnames(t): %s", colnames(t)), stderr())
    
    # The detour via a data frame below fixes a problem (issue #36) when there's only a
    # single column. Don't know why it's necessary.
    d = data.frame(t[,2:length(colnames(t))], row.names=t[,1])
    m = as.matrix(d)
    chorddiag(
      m, type = "bipartite",  groupnameFontsize =  14,
      groupColors = c(DARK_PALETTE_768X[1:length(m[,1])], LIGHT_PALETTE_768X[1:length(m[1,])]),
      categoryNames = c(
        sprintf("Organism %s", sub('^t', '', input$taxonrank)),
        sprintf(
          "Protein %s%s",
          sub('^p', '', input$proteinrank),
          ifelse( input$protstattype == 'indproteins', '', ' combinations' )
        )
      ), categorynameFontsize = 16
    )
  })
  
  output$distgraph = renderPlot({
    subc =  ifelse(
      which(PROTEIN_HIERARCHY==input$proteinrank) == length(PROTEIN_HIERARCHY),
      input$proteinrank, PROTEIN_HIERARCHY[which(PROTEIN_HIERARCHY==input$proteinrank) + 1]
    )
    d = filtered_table() %>%
      mutate(seqlen = str_length(seq)) %>%
      mutate_(
        'stat' = input$sinastat,
        'c' = input$proteinrank,
        'subc' = subc
      )
    
    ggplot(d, aes(x=c, y=stat)) + 
      geom_violin() +
      geom_sina(aes(colour=subc), method='counts') +
      scale_colour_manual(sprintf('Protein %s', sub('^p', '', subc)), values=DIV_PALETTE_768X) +
      xlab(sprintf("Protein %s", sub('^p', '', input$proteinrank))) +
      ylab('Statistic') +
      theme(
        axis.title = element_text(size=18),
        axis.text = element_text(size=12),
        axis.text.x = element_text(angle=60, hjust=1)
      )
  }) 

  output$fastaseq = downloadHandler(
    filename = 'pfitmap_sequences.faa',
    content = function(file) {
      # The session variable doesn't seem to be available in shiny-server
      #write(sprintf("DEBUG: %s: user-agent: %s", Sys.time(), session$request$HTTP_USER_AGENT), stderr())
      
      nl = "\n"
      #nl = ifelse(
        #grepl('windows', session$request$HTTP_USER_AGENT, ignore.case=T), "\r\n",
        #ifelse(grepl('mac', session$request$HTTP_USER_AGENT, ignore.case=T), "\r", "\n")
      #)
      
      d = filtered_table() %>% 
        transmute(
          taxon = sprintf("%s:%s:%s", tdomain, tphylum, tstrain),
          protein = sub('.*:', '', paste(psuperfamily, pfamily, pclass, psubclass, sep=":")),
          s = gsub('  *', '_', sprintf(">%s_%s_@%s%s%s", taxon, protein, accno, nl, seq))
        ) %>%
        select(s)
      write(d$s, file)
    },
    contentType = 'text/fasta'
  )
  
  #output$debug = renderText({
    #paste(
      #names(session),
      #names(session$request),
      #sep='\n*****\n'
    #)
  #})
  
  output$ssversion = renderText(
    (classified_proteins %>% 
      transmute(ssversion = sprintf("Source database: %s %s, downloaded %s", ss_source, ss_name, ss_version)) %>% 
      distinct())$ssversion
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
