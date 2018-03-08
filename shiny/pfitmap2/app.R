#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(feather)
library(purrr)
library(dplyr)
library(tidyr)
library(DT)
library(chorddiag)

UI_VERSION = '1.8.0'

write(sprintf("Reading feather files matching %s", Sys.getenv('PFITMAP2_FEATHER_GLOB')), stderr())
data        <- Sys.glob(Sys.getenv('PFITMAP2_FEATHER_GLOB')) %>% map(read_feather)
names(data) <- Sys.glob(Sys.getenv('PFITMAP2_FEATHER_GLOB')) %>% map(~sub('.*\\.([^.]*)\\.feather', '\\1', .))

# Defining short cut data
dbs              <- data$accessions %>% distinct(db) %>% filter(!is.na(db)) %>% arrange(db) %>% pull(db)
psuperfamilies   <- data$hmm_profiles %>% distinct(psuperfamily) %>% filter(!is.na(psuperfamily)) %>% pull(psuperfamily)
acc_taxa         <- data$taxa %>% inner_join(data$accessions, by = 'taxon')
tdomains         <- data$taxa %>% distinct(tdomain) %>% pull(tdomain)

proteins_by_db   <- data$dupfree_proteins %>%
  inner_join(data$hmm_profiles, by = 'profile') %>%
  inner_join(data$accessions, by = 'accno') %>%
  inner_join(data$taxa, by = 'taxon')

# Define UI for application 
ui <- fluidPage(
  titlePanel('pfitmap/RNRdb'),
  
  sidebarLayout(
    sidebarPanel(
      selectInput(
        'db', 'Database',
        dbs, selected = 'refseq'
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
        uiOutput('pclasses'),
        uiOutput('psubclasses')
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
      htmlOutput('ssversion'),
      textOutput('debug'),
      tabsetPanel(type = 'tabs', 
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

# Define server logic 
# For console testing:
# input <- list(psuperfamilies = 'NrdGRE', db = 'refseq', taxonrank = 'tphylum', proteinrank = 'pclass')
server <- function(input, output) {
  output$ssversion <- renderText({
    sprintf(
      "<a href='news.html'>Source database: %s %s %s, user interface version: %s</a>", 
      data$dbsources %>% pull(source), data$dbsources %>% pull(name), data$dbsources %>% pull(version), UI_VERSION
    )
  })

  # Returns a table after applying all filters the user have called for
  filtered_table = reactive({
    t <- proteins_by_db %>% filter(db == input$db)
    
    # Filters for protein hierarchy
    if ( length(input$psuperfamilies) > 0 ) t <- t %>% filter(psuperfamily %in% input$psuperfamilies)
    if ( length(input$pfamilies) > 0 )      t <- t %>% filter(pfamily %in% input$pfamilies)
    if ( length(input$pclasses) > 0 )       t <- t %>% filter(pclass %in% input$pclasses)
    if ( length(input$psubclasses) > 0 )    t <- t %>% filter(psubclass %in% input$psubclasses)

    # Filters for taxon hierarchy
    if ( length(input$tdomains) > 0 )       t <- t %>% filter(tdomain %in% input$tdomains)
    if ( length(input$tphyla) > 0 )         t <- t %>% filter(tphylum %in% input$tphyla)
    if ( length(input$tclasses) > 0 )       t <- t %>% filter(tclass %in% input$tclasses)
    if ( length(input$torders) > 0 )        t <- t %>% filter(torder %in% input$torders)
    if ( length(input$tfamilies) > 0 )      t <- t %>% filter(tfamily %in% input$tfamilies)
    if ( length(input$tgenera) > 0 )        t <- t %>% filter(tgenus %in% input$tgenera)
    if ( length(input$tspecies) > 0 )       t <- t %>% filter(tspecies %in% input$tspecies)

    return(t)
  })

  output$mainmatrix <- renderDataTable({
    ft <- filtered_table() %>% 
      group_by(!!rlang::sym(input$taxonrank), !!rlang::sym(input$proteinrank)) %>%
      summarise(n = n()) %>% ungroup() %>%
      inner_join(
        data$taxa %>% semi_join(data$accessions %>% filter(db == input$db), by = 'taxon') %>%
          group_by(!!rlang::sym(input$taxonrank)) %>%
          summarise(n_taxa = n()) %>% ungroup(),
        by = input$taxonrank
      ) %>%
      spread(!!input$proteinrank, n, fill=0)

    dt <- datatable(ft)
  })

  output$pfamilies <- renderUI({
    if ( length(input$psuperfamilies) > 0 ) {
      p <- data$hmm_profiles %>% filter(psuperfamily %in% input$psuperfamilies) %>%
        filter(!is.na(pfamily))
    } else {
      p <- data$hmm_profiles %>% distinct(pfamily) %>% 
        filter(!is.na(pfamily))
    }
    
    selectInput(
      'pfamilies', 'Families',
      p %>% distinct(pfamily) %>% pull(pfamily), multiple = T
    )
  })

  output$pclasses <- renderUI({
    if ( length(input$pfamilies) > 0 ) {
      p <- data$hmm_profiles %>% filter(pfamily %in% input$pfamilies) %>%
        filter(!is.na(pclass))
    } else {
      p <- data$hmm_profiles %>% distinct(pclass) %>% 
        filter(!is.na(pclass))
    }
    
    selectInput(
      'pclasses', 'Classes',
      p %>% distinct(pclass) %>% pull(pclass), multiple = T
    )
  })

  output$psubclasses <- renderUI({
    if ( length(input$pclasses) > 0 ) {
      p <- data$hmm_profiles %>% filter(pclass %in% input$pclasses) %>%
        filter(!is.na(psubclass))
    } else {
      p <- data$hmm_profiles %>% distinct(psubclass) %>% 
        filter(!is.na(psubclass))
    }
    
    selectInput(
      'psubclasses', 'Subclasses',
      p %>% distinct(psubclass) %>% pull(psubclass), multiple = T
    )
  })

  output$tphyla <- renderUI({
    t <- acc_taxa %>% filter(db == input$db) %>% distinct(tdomain, tphylum)
    if ( length(input$tdomains) > 0 ) t <- t %>% filter(tdomain %in% input$tdomains)
    selectInput(
      'tphyla', 'Phyla',
      t %>% arrange(tdomain, tphylum) %>% pull(tphylum),
      multiple = TRUE
    )
  })

  output$tclasses <- renderUI({
    t <- acc_taxa %>% filter(db == input$db) %>% distinct(tdomain, tphylum, tclass)
    if ( length(input$tdomains) > 0 ) t <- t %>% filter(tdomain %in% input$tdomains)
    if ( length(input$tphyla) > 0 )   t <- t %>% filter(tphylum %in% input$tphyla)
    selectInput(
      'tclasses', 'Classes',
      t %>% arrange(tdomain, tphylum, tclass) %>% pull(tclass),
      multiple = TRUE
    )
  })

  output$torders <- renderUI({
    t <- acc_taxa %>% filter(db == input$db) %>% distinct(tdomain, tphylum, tclass, torder)
    if ( length(input$tdomains) > 0 ) t <- t %>% filter(tdomain %in% input$tdomains)
    if ( length(input$tphyla) > 0 )   t <- t %>% filter(tphylum %in% input$tphyla)
    if ( length(input$tclasses) > 0 ) t <- t %>% filter(tclass %in% input$tclasses)
    selectInput(
      'torders', 'Orders',
      t %>% arrange(tdomain, tphylum, tclass, torder) %>% pull(torder),
      multiple = TRUE
    )
  })

  output$tfamilies <- renderUI({
    t <- acc_taxa %>% filter(db == input$db) %>% distinct(tdomain, tphylum, tclass, torder, tfamily)
    if ( length(input$tdomains) > 0 ) t <- t %>% filter(tdomain %in% input$tdomains)
    if ( length(input$tphyla) > 0 )   t <- t %>% filter(tphylum %in% input$tphyla)
    if ( length(input$tclasses) > 0 ) t <- t %>% filter(tclass %in% input$tclasses)
    if ( length(input$torders) > 0 )  t <- t %>% filter(torder %in% input$torders)
    selectInput(
      'tfamilies', 'Families',
      t %>% arrange(tdomain, tphylum, tclass, torder, tfamily) %>% pull(tfamily),
      multiple = TRUE
    )
  })

  output$tgenera <- renderUI({
    t <- acc_taxa %>% filter(db == input$db) %>% distinct(tdomain, tphylum, tclass, torder, tfamily, tgenus)
    if ( length(input$tdomains) > 0 )  t <- t %>% filter(tdomain %in% input$tdomains)
    if ( length(input$tphyla) > 0 )    t <- t %>% filter(tphylum %in% input$tphyla)
    if ( length(input$tclasses) > 0 )  t <- t %>% filter(tclass %in% input$tclasses)
    if ( length(input$torders) > 0 )   t <- t %>% filter(torder %in% input$torders)
    if ( length(input$tfamilies) > 0 ) t <- t %>% filter(tfamily %in% input$tfamilies)
    selectInput(
      'tgenera', 'Genera',
      t %>% arrange(tdomain, tphylum, tclass, torder, tfamily, tgenus) %>% pull(tgenus),
      multiple = TRUE
    )
  })

  output$tspecies <- renderUI({
    t <- acc_taxa %>% filter(db == input$db) %>% distinct(tdomain, tphylum, tclass, torder, tfamily, tgenus, tspecies)
    if ( length(input$tdomains) > 0 )  t <- t %>% filter(tdomain %in% input$tdomains)
    if ( length(input$tphyla) > 0 )    t <- t %>% filter(tphylum %in% input$tphyla)
    if ( length(input$tclasses) > 0 )  t <- t %>% filter(tclass %in% input$tclasses)
    if ( length(input$torders) > 0 )   t <- t %>% filter(torder %in% input$torders)
    if ( length(input$tfamilies) > 0 ) t <- t %>% filter(tfamily %in% input$tfamilies)
    if ( length(input$tgenera) > 0 )   t <- t %>% filter(tgenus %in% input$tgenera)
    selectInput(
      'tspecies', 'Species',
      t %>% arrange(tdomain, tphylum, tclass, torder, tfamily, tgenus, tspecies) %>% pull(tspecies),
      multiple = TRUE
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
