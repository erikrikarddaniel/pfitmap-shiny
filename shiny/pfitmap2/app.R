#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(pool)
library(dplyr)
library(dbplyr)

UI_VERSION = '1.8.0'

write(sprintf("Connecting to %s", Sys.getenv('PFITMAP2_SQLITE')), stderr())
dbpool <- dbPool(
  drv = RSQLite::SQLite(),
  dbname = Sys.getenv('PFITMAP2_SQLITE')
)

onStop(function() {
  poolClose(dbpool)
})

# Setting up variables to connect to some tables
accessions       <- dbpool %>% tbl('accessions')
dbsources        <- dbpool %>% tbl('dbsources')
hmm_profiles     <- dbpool %>% tbl('hmm_profiles')
dupfree_proteins <- dbpool %>% tbl('dupfree_proteins')
taxa             <- dbpool %>% tbl('taxa')

# Defining short cut data
dbs              <- accessions %>% distinct(db) %>% filter(!is.na(db)) %>% arrange(db) %>% pull(db)
psuperfamilies   <- hmm_profiles %>% distinct(psuperfamily) %>% filter(!is.na(psuperfamily)) %>% pull(psuperfamily)

proteins_by_db   <- dupfree_proteins %>%
  inner_join(hmm_profiles, by = 'profile') %>%
  inner_join(accessions, by = 'accno') %>%
  inner_join(taxa, by = 'taxon')

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
server <- function(input, output) {
  output$ssversion <- renderText({
    sprintf(
      "<a href='news.html'>Source database: %s %s %s, user interface version: %s</a>", 
      dbsources %>% pull(source),
      dbsources %>% pull(name),
      dbsources %>% pull(version),
      UI_VERSION
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

    return(t)
  })

  output$mainmatrix = renderDataTable({
    ft <- filtered_table() %>% collect()

    dt <- datatable(ft)
  })

  output$pfamilies = renderUI({
    if ( length(input$psuperfamilies) > 0 ) {
      pf <- hmm_profiles %>% filter(psuperfamily %in% input$psuperfamilies) %>%
        distinct(pfamily) %>% pull(pfamily)
    } else {
      pf <- hmm_profiles %>% distinct(pfamily) %>% pull(pfamily)
    }
    
    selectInput(
      'pfamilies', 'Families',
      pf, multiple = T
    )
  })

  output$pclasses = renderUI({
    if ( length(input$pfamilies) > 0 ) {
      pc <- hmm_profiles %>% filter(pfamily %in% input$pfamilies) %>%
        distinct(pclass) %>% pull(pclass)
    } else {
      pc <- hmm_profiles %>% distinct(pclass) %>% pull(pclass)
    }
    
    selectInput(
      'pclasses', 'Classes',
      pc, multiple = T
    )
  })

  output$psubclasses = renderUI({
    if ( length(input$pclasses) > 0 ) {
      pc <- hmm_profiles %>% filter(pclass %in% input$pclasses) %>%
        distinct(psubclass) %>% pull(psubclass)
    } else {
      pc <- hmm_profiles %>% distinct(psubclass) %>% pull(psubclass)
    }
    
    selectInput(
      'psubclasses', 'Subclasses',
      pc, multiple = T
    )
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
