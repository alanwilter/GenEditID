shinyServer(function(input, output, session){
  
  # selected entities
  
  selection <- reactiveValues(
    selectedProject = NULL
  )
  
  
#### PROJECT TAB
  #Project table
  output$ProjectTable <- DT::renderDataTable(
    datatable(as.data.frame(project),
              rownames = FALSE,
              selection = "single",
              extensions = "Buttons",
              options = list(
                dom = 'Bfrtip',
                buttons = I('pageLength'),
                pageLength = 10,
                lengthMenu = list(c(5, 10, 25, -1), c('Show 5 rows', 'Show 10 rows', 'Show 25 rows', 'Show all rows')),
                searchHighlight = TRUE
              )
    ),
    server = TRUE
  )
  
  #Create proxy object to manipulate project table
  ProjectTableProxy <- dataTableProxy('ProjectTable')
  
  observeEvent(selection$selectedProject, {
    selectRows(
      ProjectTableProxy,
    which(as.data.frame(project)$project_geid == selection$selectedProject)
    )
  })
  
  observe({
    selectedRow <- input$projectTable_rows_selected
    if (!is.null(selectedRow))
    {
      selectedProject <- project %>%
        slice(selectedRow) %>%
        select(project_geid) %>%
        unlist %>%
        as.character
      selection$selectedProject <- selectedProject
    }
  })
  
  #Guide table
  output$GuideTable <- renderTable(as.data.frame(guide)[, -c(1:2)])
  
  
#### EXPLORATORY ANALYSIS TAB
  output$proteinPlot <- renderPlotly(protein_abundance_plot)
  output$growthcurvePlot <- renderPlotly(clone_growth_curve)
  
#### RESULTS TAB
  
  
  
})

            














