library(dplyr)
library(ggplot2)
library(svglite)
library(RPostgreSQL)

fun.loadPlateData <- function(conn, layoutId) {
  sql <- paste0(
    'select cl.name as cell_line_name, c.name as clone_name, t.name as target_name,
        g.name as guide_name, wc.content_type, w.row, w.column, vr.ge_score
        from well w
        inner join experiment_layout el on w.experiment_layout_id = el.id
        left join well_content wc on w.well_content_id = wc.id
        left join sequencing_library_content slc on slc.well_id = w.id
        left join variant_result vr on vr.sequencing_library_content_id = slc.id
        left join clone c on wc.clone_id = c.id
        left join cell_line cl on c.cell_line_id = cl.id
        left join guide_well_content_association gwca on gwca.well_content_id = wc.id
        inner join guide g on gwca.guide_id = g.id
        inner join target t on g.target_id = t.id
        where wc.content_type is not null
        and el.geid = \'', layoutId, '\'')
  
  cat(sql)
  cat('\n')
  
  query <- as.data.frame(dbGetQuery(conn, sql))
  
  query$classifier <- create_classifier(query$cell_line_name, query$clone_name, query$guide_name, query$content_type)
  
  print(query)
  
  query
}

driver <- dbDriver("PostgreSQL")
conn <- dbConnect(driver, user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")

fun.layoutplot <- function(conn, layoutId){
    
    tryCatch(
      {
        plateData <- fun.loadPlateData(conn, layoutId)
        layout_plot <- plateData %>%
          as.data.frame %>%
          mutate(ge_score = coalesce(ge_score, 0L)) %>% # replace NA's with 0's
          plot_ly(x = ~column, y = ~row,
                  mode = "markers", type = "scatter", color = ~ge_score, colors = c('grey90', 'orange', 'green4'),
                  marker = list(size = 15),
                  # Hover text:
                  text = ~paste0(row, column, ' score: ', ge_score)
          ) %>%
          layout(yaxis = list(autorange = "reversed"))
      return(layout_plot)
        
      },
      finally =
      {
        dbDisconnect(conn)
      }
    )
}