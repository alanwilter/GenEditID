shinyUI(fluidPage(
  titlePanel(h1(strong('ICW and IncuCyte'), style='color:steelblue')),
  br(),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel('Controls',
                 p(h4('Primer design:', style='color:steelblue')),
                 
                 
                 
                 fileInput("fileinput", "Choose CSV File",
                            # accept = c(
                             #  "text/csv",
                              # "text/comma-separated-values,text/plain",
                               #".csv"),
                             multiple = T
                   
                 )
                )
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel('Plots',
                  fluidRow(
                    column(1, br(), br(), br(), downloadButton('plotICWDownload', ''), br(), br()
                    ),
                    column(11, h3('IncellWestern'),
                           plotlyOutput('ICWplot'))
                  ),
                 fluidRow(
                   column(1, br(), br(), br(), downloadButton('plotInRDownload', ''), br(), br()),
                   column(11, h3('IncuCyte'),
                          plotlyOutput('Incucyteplot')),
                          tableOutput('tablecsv')
                 )
                  ,
                  fluidRow(
                    column(1, br(), br(), br(), downloadButton('plotInCDownload', ''), br(), br()),
                    column(11, h3('IncuCyte'),
                           plotlyOutput('IncucyteCurveplot')) #,
                    #tableOutput('tablecsv')
                  )
        )
      )
    )
  )
  
))