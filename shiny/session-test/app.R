#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

ui <- fluidPage(
  titlePanel("Testing the session variable"),
   
  flowLayout(
   textOutput('textarea')
  )
)

server <- function(input, output, session) {
  output$textarea = renderText({
    #write(sprintf("DEBUG: %s: user-agent: %s", Sys.time(), session$request$HTTP_USER_AGENT), stderr())
    sprintf("HTTP_USER_AGENT: %s", session$request$HTTP_USER_AGENT)
  })
}

shinyApp(ui = ui, server = server)