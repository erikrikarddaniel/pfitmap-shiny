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

# Reading data and transforming
write(sprintf("LOG: %s: Reading data files", Sys.time()), stderr())

con <- DBI::dbConnect(RSQLite::SQLite(), '../../hmmer-import/R-test/hmmsearch2classification.03.sqlite3')
bestscoring <- con %>% tbl('bestscoring')
hmm_profiles <- con %>% tbl('hmm_profiles')

write(sprintf("LOG: %s: Connect done", Sys.time()), stderr())

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel('TEST'),
  
  sidebarLayout(
    sidebarPanel('sidebar'),
    mainPanel(
      'main',
      plotOutput('classhistogram')
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$classhistogram <- renderPlot({
    hmm_profiles %>%
      inner_join(bestscoring, by = 'profile') %>% 
      group_by(psubclass) %>%
      summarise(n = n()) %>%
      ungroup() %>%
      collect() %>%
      ggplot(aes(x = psubclass, y = n)) +
      geom_col()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)