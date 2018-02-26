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
  poolClose(pool)
})

dbsources <- dbpool %>% tbl('dbsources')

# Define UI for application 
ui <- fluidPage(
  titlePanel('pfitmap/RNRdb'),
  
  sidebarLayout(
    sidebarPanel(
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
}

# Run the application 
shinyApp(ui = ui, server = server)
