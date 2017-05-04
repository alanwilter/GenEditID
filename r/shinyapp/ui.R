shinyUI(fluidPage(
  titlePanel(h1(strong('GE Project Analysis'), style='color:steelblue')),
  br(),
  
  tabsetPanel(
    
    
    
    tabPanel('Project',
      p(h4('Project selection:', style='color:steelblue')),
      DT::dataTableOutput('ProjectTable'),
      br(),
      p('sgRNA guide details:', style='color:steelblue'),
      tableOutput('GuideTable'),
      br(),
      p('Project progress:', style='color:steelblue'),
      p('show Growth, protein and NGS slots, reacting to presence/absence of data with colours')
    ),
  
    
    
    tabPanel('Exploratory Analysis',
      p(),
      sidebarLayout(
        
        sidebarPanel(width = 2,
            radioButtons("EAradio", label = h3("Select dataset"),
                choices = list("Growth" = 1, "Protein" = 2, "NGS" = 3, 'Combined data' = 4), 
                selected = 1)            
        ),
        mainPanel(
            conditionalPanel(
               condition = "input.EAradio == 1",
               column(6, plotlyOutput('growthcurvePlot', height = "600px")),
               column(6, plotlyOutput('growthratePlot', height = "800px"))
            ),
            conditionalPanel(
               condition = "input.EAradio == 2",
               fillPage(plotlyOutput('proteinPlot', height = "800px"))
            ),
            # conditionalPanel(
            #    condition = "input.EAradio == 3",
            #    plotlyOutput()
            # ),
            conditionalPanel(
               condition = "input.EAradio == 4",
               fillPage(plotlyOutput('combinedPlot', height = "800px"))
            )
                  
                  
                  # fluidRow(
                  #   column(1, br(), br(), br(), downloadButton('plotICWDownload', ''), br(), br()
                  #   ),
                  #   column(11, h3('IncellWestern'),
                  #          plotlyOutput('ICWplot'))
                  # )
        )
      )
    ),
    
    
    
    tabPanel('Results',
      p(),
      sidebarLayout(
        
        sidebarPanel('Clones'),
        
        mainPanel('Plots')
      )
    )
    
  )
  

))
  
  
#   sidebarLayout(
#     sidebarPanel(
#   
#         
#         
#                  
#                 )
#       )
#     ),
#     mainPanel(
#       tabsetPanel(
#         
#         tabPanel('Exploratory plots'),
#         
#         tabPanel('Results'),
#         
#         tabPanel('Plots',
#                   ,
#                  fluidRow(
#                    column(1, br(), br(), br(), downloadButton('plotInRDownload', ''), br(), br()),
#                    column(11, h3('IncuCyte'),
#                           plotlyOutput('Incucyteplot')),
#                           tableOutput('tablecsv')
#                  )
#                   ,
#                   fluidRow(
#                     column(1, br(), br(), br(), downloadButton('plotInCDownload', ''), br(), br()),
#                     column(11, h3('IncuCyte'),
#                            plotlyOutput('IncucyteCurveplot')) #,
#                     #tableOutput('tablecsv')
#                   )
#         )
#       )
#     )
#   )
#   
# ))