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

# Defining short cut data
dbs            <- accessions %>% distinct(db) %>% filter(!is.na(db)) %>% arrange(db) %>% pull(db)
psuperfamilies <- hmm_profiles %>% distinct(psuperfamily) %>% filter(!is.na(psuperfamily)) %>% pull(psuperfamily)

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
      htmlOutput('ssversion')
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
