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

# Read all feather files into a common list object.
write(paste0(Sys.time(), ": Reading feather files matching ", Sys.getenv('PFITMAP2_FEATHER_GLOB')), stderr())
data        <- Sys.glob(Sys.getenv('PFITMAP2_FEATHER_GLOB')) %>% map(read_feather)
names(data) <- Sys.glob(Sys.getenv('PFITMAP2_FEATHER_GLOB')) %>% map(~sub('.*\\.([^.]*)\\.feather', '\\1', .))

# Defining shortcut data for easy access to the various parts of the data object
write(paste0(Sys.time(), ": Done reading data, setting up global defines"), stderr())
dbs              <- data$accessions %>% distinct(db) %>% filter(!is.na(db)) %>% arrange(db) %>% pull(db)
psuperfamilies   <- data$hmm_profiles %>% distinct(psuperfamily) %>% filter(!is.na(psuperfamily)) %>% pull(psuperfamily)
acc_taxa         <- data$taxa %>% inner_join(data$accessions, by = 'taxon')
tdomains         <- data$taxa %>% distinct(tdomain) %>% pull(tdomain)

proteins_by_db   <- data$dupfree_proteins %>%
  inner_join(data$hmm_profiles, by = 'profile') %>%
  inner_join(data$accessions, by = 'accno') %>%
  inner_join(data$taxa, by = 'taxon')

write(paste0(Sys.time(), ": Done defining data"), stderr())

# Define UI for application 
ui <- fluidPage(
  titlePanel('pfitmap/RNRdb'),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        'db', 'Database',
        dbs, selected = 'refseq'
      )
    ),
    mainPanel(
      htmlOutput('ssversion')
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
}

# Run the application 
shinyApp(ui = ui, server = server)
