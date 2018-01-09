#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggforce)
library(stringr)
library(DT)
library(chorddiag)

# Some constants
PROFILES_VERSION = '0.7'
UI_VERSION = '1.1.0'

PROTEIN_HIERARCHY = c( 'psuperfamily', 'pfamily', 'pclass', 'psubclass', 'pgroup' )
TAXON_HIERARCHY = c( 'tdomain', 'tkingdom', 'tphylum', 'tclass', 'torder', 'tfamily', 'tgenus', 'tspecies', 'tstrain' )

TRAIT_PRESENCE_BOTH = 'Present/absent'
TRAIT_PRESENCE_ONLY_PRESENT = 'Only presences'
TRAIT_PRESENCE_ONLY_ABSENT = 'Only absences'

INDPROTEINS  = 'indproteins'
COMBPROTEINS = 'combproteins'
COMBDOMAINS  = 'combdomains'

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
write(sprintf("LOG: %s: Opening the database %s", Sys.time(), Sys.getenv('PFITMAP_SQLITEDB')), stderr())
db <- DBI::dbConnect(RSQLite::SQLite(), Sys.getenv('PFITMAP_SQLITEDB'))
write(sprintf("LOG: %s: Database opened", Sys.time()), stderr())

# Define UI
ui <- fluidPage(
  titlePanel('pfitmap/RNRdb'),
  
  sidebarLayout(
    sidebarPanel( 
      radioButtons(
        'protstattype', 'Type of protein statistic',
        list(
          'Individual proteins' = INDPROTEINS,
          'Combinations of proteins' = COMBPROTEINS,
          'Combinations of domains' = COMBDOMAINS
        ),
        selected = INDPROTEINS
###      ),
###      selectInput(
###        'db', 'Database',
###        dbs, selected = 'ref'
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
###      ),
###      wellPanel(
###        selectInput(
###          'psuperfamilies', 'Protein superfamilies',
###          c('', psuperfamilies), selected = c('NrdGRE'), 
###          multiple=T
###        ),
###        uiOutput('pfamilies'),
###        uiOutput('pclasses')
###      ),
###      wellPanel(
###        selectInput(
###          'tdomains', 'Taxonomic domains',
###          c('', tdomains), multiple = T
###        ),
###        uiOutput('tphyla'),
###        uiOutput('tclasses'),
###        uiOutput('torders'),
###        uiOutput('tfamilies'),
###        uiOutput('tgenera'),
###        uiOutput('tspecies')
      )
    ),
    mainPanel(
      htmlOutput('ssversion'),
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
          downloadLink('fastaseq', 'Download sequences in fasta format'),
          htmlOutput('sequencelist')
        )
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # Reactive methods outside assignments

  # Filters the taxa table based on all selectors. Don't use to set
  # lists for input boxes; works for the select in filtered_table().
  filtered_taxa = reactive({
    t = taxa %>% filter(db == input$db)

    if ( length(input$tdomains) > 0 )  t = t %>% filter(tdomain  %in% input$tdomains)
    if ( length(input$tphyla) > 0 )    t = t %>% filter(tphylum  %in% input$tphyla)
    if ( length(input$tclasses) > 0 )  t = t %>% filter(tclass   %in% input$tclasses)
    if ( length(input$torders) > 0 )   t = t %>% filter(torder   %in% input$torders)
    if ( length(input$tfamilies) > 0 ) t = t %>% filter(tfamily  %in% input$tfamilies)
    if ( length(input$tgenera) > 0 )   t = t %>% filter(tgenus   %in% input$tgenera)
    if ( length(input$tspecies) > 0 )  t = t %>% filter(tspecies %in% input$tspecies)

    t
  })

  # Filters proteins similar to the above.
  filtered_proteins = reactive({
    p = classified_proteins %>% filter(db == input$db)

    if ( length(input$psuperfamilies) > 0 )  p = p %>% filter(psuperfamily  %in% input$psuperfamilies)
    if ( length(input$pfamilies) > 0 )       p = p %>% filter(pfamily       %in% input$pfamilies)
    if ( length(input$pclasses) > 0 )        p = p %>% filter(pclass        %in% input$pclasses)

    p
  })
  
  # Returns a table after applying all filters the user have called for
  filtered_table = reactive({
    ### WORK IN PROGRESS
    t = db %>% tbl('accessions') %>% 
      filter(db == input$db) %>%
      inner_join(db %>% tbl('bestscoring'), by = 'accno') %>%
      inner_join(db %>% tbl('hmm_profiles'), by = 'profile') %>%
      collect() %>%                    # Database does not support rank()
      group_by(accno, taxon) %>% mutate(r = rank(accto)) %>% ungroup() %>% filter(r == 1) %>% # One count per db, taxon and real sequence
      mutate(accno = accto) %>% select(-accno, -r) %>%   # accto is confusing comes from tblout processing
      inner_join(filtered_taxa() %>% select(db, ncbi_taxon_id), by=c('db', 'ncbi_taxon_id'))
    
    # Filters for protein hierarchy
    if ( length(input$psuperfamilies) > 0 ) { t = t %>% filter(psuperfamily %in% input$psuperfamilies) }
    if ( length(input$pfamilies) > 0 ) { t = t %>% filter(pfamily %in% input$pfamilies) }
    if ( length(input$pclasses) > 0 ) { t = t %>% filter(pclass %in% input$pclasses) }

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
    
    ###write(sprintf("DEBUG: %s: rows in table: %d", Sys.time(), length(t[,1])), stderr())
    
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
    
    ###write(sprintf("DEBUG: %s: rows in grouped table: %d", Sys.time(), length(d[,1])), stderr())
    
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

  # Calculates all present combinations of proteins at the proteinrank selected,
  # groups by the combinations, combines with n_genomes from taxa and returns.
  combdomains_sums_table = reactive({
    ###write(sprintf("combdomains_sums_table, protstattype: %s", input$protstattype), stderr())
    d = filtered_table() %>%
      select_('accno', 'tsort', 'tcolour', 'taxon_tooltip', 'ncbi_taxon_id', input$taxonrank, input$proteinrank) %>%
      distinct()
    names(d)[names(d)==input$proteinrank] = 'domain'
    d = d %>%
      union(
        domain_hits %>% select(accno, domain) %>%  
          inner_join(d %>% select(-domain), by='accno')
      ) %>%
      group_by_('tsort', 'tcolour', 'taxon_tooltip', 'accno', input$taxonrank) %>%
      mutate(n = row_number(domain)) %>%
      ungroup() %>%
      unite(domain, n, domain) %>%
      spread(domain, domain, fill='')
    d = d %>% unite(comb, 7:length(colnames(d0)), sep=':') %>%
      mutate(
        comb = gsub(
          '[0-9][0-9]*_', '', gsub(
            '::*', '', gsub(
              '::*$', '', gsub(
                '^::*', '', comb
              )
            )
          )
        )
      ) %>%
      group_by_('tsort', 'tcolour', 'taxon_tooltip', 'comb', input$taxonrank) %>%
      summarise(n = n()) %>%
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
    pf = classified_proteins %>% filter(db == input$db)

    if ( length(input$psuperfamilies) > 0 ) pf = pf %>% filter(psuperfamily %in% input$psuperfamilies)

    pf = (pf %>% select(pfamily) %>% distinct())$pfamily
    
    selectInput(
      'pfamilies', 'Families',
      pf, multiple = T
    )
  })
  
  output$pclasses = renderUI({
    pc = classified_proteins %>% filter(db == input$db)

    if ( length(input$psuperfamilies) > 0 ) pc = pc %>% filter(psuperfamily %in% input$psuperfamilies)
    if ( length(input$pfamilies) > 0) pc = pc %>% filter(pfamily %in% input$pfamilies)

    pc = (pc %>% select(pclass) %>% distinct() %>% arrange(pclass))$pclass
    
    selectInput(
      'pclasses', 'Classes',
      pc, multiple = T
    )
  })

  output$tphyla = renderUI({
    ###write(sprintf("input$tdomains, len: %d: %s", length(input$tdomains), input$tdomains), stderr())
    # Calling the filtered_taxa() function doesn't work. It seems the
    # reference to lower taxonomic rank inputs resets everything.
    t = taxa %>% filter(db == input$db)
    if ( length(input$tdomains) > 0 )  t = t %>% filter(tdomain  %in% input$tdomains)
    t = t %>% inner_join(filtered_proteins() %>% select(ncbi_taxon_id), by='ncbi_taxon_id')
    selectInput(
      'tphyla', 'Phyla',
      (t %>% select(tphylum) %>% distinct() %>% arrange(tphylum))$tphylum,
      multiple = T
    )
  })
  
  output$tclasses = renderUI({
    ###write(sprintf("tphyla: %d: %s", length(input$tphyla), input$tphyla), stderr())
    t = taxa %>% filter(db == input$db)

    if ( length(input$tdomains) > 0 )  t = t %>% filter(tdomain  %in% input$tdomains)
    if ( length(input$tphyla) > 0 )    t = t %>% filter(tphylum  %in% input$tphyla)

    t = t %>% inner_join(filtered_proteins() %>% select(ncbi_taxon_id), by='ncbi_taxon_id')

    selectInput(
      'tclasses', 'Classes',
      (t %>% select(tclass) %>% distinct() %>% arrange(tclass))$tclass,
      multiple = T
    )
  })
  
  output$torders = renderUI({
    ###write(sprintf("tclasses: %d: %s", length(input$tclasses), input$tclasses), stderr())
    t = taxa %>% filter(db == input$db)
    if ( length(input$tdomains) > 0 )  t = t %>% filter(tdomain  %in% input$tdomains)
    if ( length(input$tphyla) > 0 )    t = t %>% filter(tphylum  %in% input$tphyla)
    if ( length(input$tclasses) > 0 )  t = t %>% filter(tclass   %in% input$tclasses)
    t = t %>% inner_join(filtered_proteins() %>% select(ncbi_taxon_id), by='ncbi_taxon_id')
    selectInput(
      'torders', 'Orders',
      (t %>% select(torder) %>% distinct() %>% arrange(torder))$torder,
      multiple = T
    )
  })
  
  output$tfamilies = renderUI({
    ###write(sprintf("torders: %d: %s", length(input$torders), input$torders), stderr())
    t = taxa %>% filter(db == input$db)
    if ( length(input$tdomains) > 0 )  t = t %>% filter(tdomain  %in% input$tdomains)
    if ( length(input$tphyla) > 0 )    t = t %>% filter(tphylum  %in% input$tphyla)
    if ( length(input$tclasses) > 0 )  t = t %>% filter(tclass   %in% input$tclasses)
    if ( length(input$torders) > 0 )   t = t %>% filter(torder   %in% input$torders)
    t = t %>% inner_join(filtered_proteins() %>% select(ncbi_taxon_id), by='ncbi_taxon_id')
    selectInput(
      'tfamilies', 'Families',
      (t %>% select(tfamily) %>% distinct() %>% arrange(tfamily))$tfamily,
      multiple = T
    )
  })
  
  output$tgenera = renderUI({
    ###write(sprintf("tfamilies: %d: %s", length(input$tfamilies), input$tfamilies), stderr())
    t = taxa %>% filter(db == input$db)
    if ( length(input$tdomains) > 0 )  t = t %>% filter(tdomain  %in% input$tdomains)
    if ( length(input$tphyla) > 0 )    t = t %>% filter(tphylum  %in% input$tphyla)
    if ( length(input$tclasses) > 0 )  t = t %>% filter(tclass   %in% input$tclasses)
    if ( length(input$torders) > 0 )   t = t %>% filter(torder   %in% input$torders)
    if ( length(input$tfamilies) > 0 ) t = t %>% filter(tfamily  %in% input$tfamilies)
    t = t %>% inner_join(filtered_proteins() %>% select(ncbi_taxon_id), by='ncbi_taxon_id')
    selectInput(
      'tgenera', 'Genera',
      (t %>% select(tgenus) %>% distinct() %>% arrange(tgenus))$tgenus,
      multiple = T
    )
  })
  
  output$tspecies = renderUI({
    t = taxa %>% filter(db == input$db)
    if ( length(input$tdomains) > 0 )  t = t %>% filter(tdomain  %in% input$tdomains)
    if ( length(input$tphyla) > 0 )    t = t %>% filter(tphylum  %in% input$tphyla)
    if ( length(input$tclasses) > 0 )  t = t %>% filter(tclass   %in% input$tclasses)
    if ( length(input$torders) > 0 )   t = t %>% filter(torder   %in% input$torders)
    if ( length(input$tfamilies) > 0 ) t = t %>% filter(tfamily  %in% input$tfamilies)
    if ( length(input$tgenera) > 0 )   t = t %>% filter(tgenus   %in% input$tgenera)
    t = t %>% inner_join(filtered_proteins() %>% select(ncbi_taxon_id), by='ncbi_taxon_id')
    selectInput(
      'tspecies', 'Species',
      (t %>% select(tspecies) %>% distinct() %>% arrange(tspecies))$tspecies,
      multiple = T
    )
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
        combproteins = combproteins_sums_table() %>% spread(comb, n, fill=0),
        combdomains  = combdomains_sums_table() %>% spread(comb, n, fill=0)
      )
      ###write(sprintf("DEBUG: %s: \tcolnames: %s", Sys.time(), paste(colnames(t), collapse=", ")), stderr())
      #write.csv(t, stderr(), row.names=F)
      
      if ( input$taxonomysort ) {
        t = t %>% arrange(tsort)
      }
      
      # This is to get the right column names, a bit involved perhaps...
      t = t %>% mutate_('Taxon'=input$taxonrank, `N. genomes`='n_genomes') %>%
        mutate(Taxon = sprintf("<span title='%s'>%s</span>", taxon_tooltip, Taxon))
      cnames = colnames(t)
      ###write(sprintf("DEBUG: %s: \tcolnames: %s", Sys.time(), paste(colnames(t), collapse=", ")), stderr())
      
      # Colours for heatmap
      brks = quantile(c(0,1), probs = seq(.0, 1, .05), na.rm = TRUE)
      clrs = round(seq(255, 40, length.out = length(brks) + 1), 0) %>% { paste0("rgb(255,", ., ",", ., ")") }
      
      # Hide some columns
      invisible = c(c(0, 1, 2, 3), grep('fraction', cnames) - 0)
      ###write(sprintf("--> invisible: %s", paste(invisible, collapse=", ")), stderr())
      ###write(sprintf("--> invisible names: %s", paste(cnames[invisible], collapse=", ")), stderr())

      ###write(sprintf("DEBUG: %s: \tcolnames: %s", Sys.time(), paste(colnames(t), collapse=", ")), stderr())
      t = t %>%
        select(tcolour, c(length(cnames)-1,length(cnames),8:length(cnames)-2))
      ###write(sprintf("DEBUG: %s: \tcolnames: %s", Sys.time(), paste(colnames(t), collapse=", ")), stderr())
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
      
      for (ff in grep('fraction', cnames, value = TRUE)) {
        dt = dt %>% formatStyle(
          sub('fraction', '', ff), ff,
          backgroundColor = styleInterval(brks, clrs)
        )
      }
      ###write(sprintf("DEBUG: %s: \trenderDataTable done", Sys.time()), stderr())
      
      dt
    }
  )
  
  output$chordgraph = renderChorddiag({
    t = switch(
      input$protstattype,
      indproteins  = indproteins_sums_table() %>% spread(proteinrank, n, fill=0),
      combproteins = combproteins_sums_table() %>% spread(comb, n, fill=0),
      combdomains  = combdomains_sums_table() %>% spread(comb, n, fill=0)
    ) %>% 
      select(-tsort, -tcolour, -taxon_tooltip, -n_genomes)
    
    # Delete all columns ending in 'fraction'
    for ( c in colnames(t)[grep('fraction', colnames(t))] ) {
      t = t %>% select_(sprintf("-`%s`", c))
    }
    
    ###write_tsv(t, '/tmp/chorddiag.t.tsv.gz')
    ###write(sprintf("colnames(t): %s", paste(colnames(t), collapse = ', ')), stderr())
    ###write(sprintf("row names: %s", paste(t[,1], collapse = ', ')), stderr())
    
    # The detour via a data frame below fixes a problem (issue #36) when there's only a
    # single column. Don't know why it's necessary.
    d = data.frame(t[, sapply(t, is.numeric)])
    rownames(d) = t %>% pull(tsort)
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
          accno, nl, seq
        ) %>%
        group_by(taxon, protein, seq) %>%
        summarise(accnos = paste(accno, collapse="_")) %>%
        ungroup() %>%
        transmute(
          s = gsub('  *', '_', sprintf(">%s_%s_@%s%s%s", taxon, protein, accnos, nl, seq))
        )
      write(d$s, file)
    },
    contentType = 'text/fasta'
  )

  output$sequencelist = renderText({
    d = filtered_table() %>%
      transmute(
        db, ncbi_taxon_id, tstrain,
        acclink = ifelse(
          db == 'pdb',
          sprintf("<a href='http://www.rcsb.org/pdb/explore/explore.do?structureId=%s'>%s</a>", sub('_.*', '', accno), accno),
          ifelse(
            db %in% c('ref', 'gb'),
            sprintf("<a href='https://www.ncbi.nlm.nih.gov/protein/%s'>%s</a>", accno, accno),
            ifelse(
              db == 'sp', 
              sprintf("<a href='http://www.uniprot.org/uniprot/%s'>%s</a>", accno, accno),
              sprintf("%s:%s", db, accno)
            )
          )
        ),
        taxon = sprintf(
          "<a href='https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=%s&lvl=3&lin=f&keep=1&srchmode=1&unlock'>%s:%s:%s</a>", 
          ncbi_taxon_id, tdomain, tphylum, tstrain
        ),
        protein = sub('.*:', '', paste(psuperfamily, pfamily, pclass, psubclass, sep=":")),
        seq
      ) %>%
      group_by(db, ncbi_taxon_id, tstrain, taxon, protein, seq) %>%
      summarise( acclinks = paste(acclink, collapse=", ")) %>%
      ungroup() %>%
      arrange(tstrain, protein) %>%
      transmute(s = sprintf(">%s %s (%s)\n%s", taxon, protein, acclinks, seq))

    sprintf("<pre>%s</pre>", paste(d$s, sep="\n", collapse="\n"))
  })
  
  #output$debug = renderText({
    #paste(
      #names(session),
      #names(session$request),
      #sep='\n*****\n'
    #)
  #})
  
  output$ssversion = renderText({
    sprintf(
      "<a href='news.html'>%s</a>",
      sv = paste(
        sprintf("Source database: %s.", db %>% tbl('dbsources') %>% collect() %>% unite(s, 1:3, sep = ' ')),
        sprintf("Profiles version %s.", PROFILES_VERSION),
        sprintf("UI version %s.", UI_VERSION)
      )
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
