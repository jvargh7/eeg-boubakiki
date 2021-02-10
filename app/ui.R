
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("EEG Analysis of Bouba Kiki Experiment- SJRI"),
  sidebarLayout(
    sidebarPanel(
      fileInput('edf_file', 'Choose EDF file',
                accept=c('edf', 'EEG Data Files', '.edf')),
      sliderInput("averaging","Select milliseconds to average for baseline correction", min = 0,max=1000,value=400,ticks=TRUE,width='400px'),
      sliderInput("stim_start","Select milliseconds to average for baseline correction", min = 0,max=1000,value=150,ticks=TRUE,width='400px'),
      sliderInput("stim_stop","Select milliseconds to average for baseline correction", min = 0,max=2000,value=1000,ticks=TRUE,width='400px'),
      actionButton("goButton", "Go!"),
      
      downloadButton('downloadresult',label='Download Result')
    ),
    mainPanel(
      textOutput("statement"),
      selectInput("channel", label = "Select Channel",choices = c("F3","Fz","F4","C3","Cz","C4","P3","P4")),
      sliderInput("range","Select Range to check for N400", min = -150,max=1000,value=c(350,550),ticks=TRUE,width='400px',dragRange=TRUE),
      fluidRow(
        column(6,tableOutput('result_table'))
      ),
      plotOutput("channel_plot")
       
    )
    
  )

))


