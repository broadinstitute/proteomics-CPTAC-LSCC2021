#################################################################
## Filename: ui.R
## Created: March 13, 2020
## Author(s): Karsten Krug
## Purpose: Shiny-app to visualize data from CPTAC LSCC discovery cohort
## This file defines the user interface of the Shiny-app.
#################################################################
library(shiny)

########################################################
## Define UI
########################################################
shinyUI(fluidPage(
    
    ## Application title
    titlePanel(HTML(TITLESTRING), windowTitle = WINDOWTITLE),
    
    HTML('<br>'),
    HTML('<b>Supplemental data:</b>'),
    
    fluidRow(
      
      ##################################################
      ## LEFT panel
      ##################################################
      column(3, wellPanel(
        
        ## authentification
        uiOutput("auth.user"),
        
        ## input
        uiOutput("ui.input")
      )
      
      ),
      
      ##################################################
      ## RIGHT panel
      ##################################################
      column(9,
             plotOutput("plot")
      ) ## end column
      
    ) ## end fluiRow
))
