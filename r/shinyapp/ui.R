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
               column(6, plotlyOutput('growthcurvePlot')),
               column(6, plotlyOutput('growthratePlot'))
            ),
            conditionalPanel(
               condition = "input.EAradio == 2",
               plotlyOutput('proteinPlot', height = "600px")
            ),
            conditionalPanel(
               condition = "input.EAradio == 3",
               tabsetPanel(
                 tabPanel('Indel positional ranges',
                          plotlyOutput('NGS.indelrangesPlot')
                 ),
                 tabPanel('Proportions',
                    column(5, plotlyOutput('NGS.mutationsPlot', height = '200px'),
                              plotlyOutput('NGS.distancetocutsitePlot', height = '200px'),
                              plotlyOutput('NGS.zygosityPlot', height = '200px')
                           ),
                    column(7, plotlyOutput('NGS.variantsPlot', height = '200px'),
                              plotlyOutput('NGS.allelesPlot')
                           )
                 )
               )
            ),
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
        mainPanel('Plots',
                  plotlyOutput('plate.scoresPlot', height = '200px', width = '400px')
        )
      )
    )
    
  )
))
