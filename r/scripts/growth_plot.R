# DESCRIPTION:
  # R script to process growth data from the IncuCyte.
  # Functions to read from database, calculate growth slopes using the grofit package and produce plots of growth curves and growth rates.

library(dplyr)
library(grofit)
library(plotly)

fun.growth_readDB <- function() {
# connect to database
    db <- src_postgres(user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")
    
    # connect to tables
    clone <- tbl(db, "clone") %>% rename(clone_id=id, clone_name=name)
    well <- tbl(db, "well") %>% rename(well_id=id)
    well_content <- tbl(db, "well_content") %>% rename(well_content_id=id)
    plate <- tbl(db, "plate") %>% rename(plate_id=id)
    growth <- tbl(db, "growth") %>% rename(growth_id=id)
    target <- tbl(db, "target") %>% rename(target_id=id, target_name=name)
    guide <- tbl(db, "guide") %>% rename(guide_id=id, guide_name=name)
    guide_well_content_association <- tbl(db, "guide_well_content_association")
    
    # read clone growth data
    clone_growth_data <- left_join(well, growth, by="well_id") %>%
      inner_join(., well_content, by="well_content_id") %>%
      inner_join(., clone, by="clone_id") %>%
      inner_join(., cell_line, by="cell_line_id") %>%
      inner_join(., plate, by="plate_id") %>%
      left_join(., guide_well_content_association, by="well_content_id") %>%
      left_join(., guide, by="guide_id") %>%
      left_join(., target, by="target_id") %>%
      collect %>%
      mutate(position=paste0(row, column)) %>%
      mutate(classifier=create_classifier(cell_line_name, clone_name, guide_name, content_type)) %>%
      rename(Content = classifier, Plate=plate_id, Target=target_name, Guide=guide_name, Well=position, Elapsed=hours, Confluence=confluence_percentage) %>%
      filter(!is.na(Confluence)) %>%
      mutate(Plate = as.factor(Plate), Well = as.factor(Well), Guide = as.factor(Guide))
    return(clone_growth_data)
}

#-----------------------------
# CALCULATIONS

# grofit -Process clone_growth_data to use with grofit
fun.growth_calcslope <- function(clone_growth_data) {
  
    # Calculate growth increase between time 93 and time 0.
    growthincrease <- clone_growth_data[clone_growth_data$Elapsed == 93, 'Confluence'] - clone_growth_data[clone_growth_data$Elapsed == 0, 'Confluence']
    colnames(growthincrease) <- 'Growth.increase'
    growthincrease <- data.frame(clone_growth_data[clone_growth_data$Elapsed == 0, c('Plate', 'Well')], growthincrease)
    
    # Trim out wells with growth increase (time 93 - time 0) < 12%
      # it could be worth it to make the growthincrease threshold manually selectable (so a user can say 'I prefer faster growers, etc')
    growthincrease <- growthincrease[growthincrease$Growth.increase > 12,] %>%
      mutate(index = paste0(Plate, Well))
    
    # create dataset (data.frame) to use with grofit (grofit requires 3 columns before the numeric data, otherwise it will fail)
    clone_growth_data_rate <- select(clone_growth_data, Plate, Well, Content, Elapsed, Confluence) %>%
      mutate(classifier = as.factor(Content)) %>%
      filter(paste0(clone_growth_data$Plate, clone_growth_data$Well) %in% growthincrease$index) %>%
      dcast(Plate + Well + Content ~ Elapsed, value.var = 'Confluence')
    
    # create timeset (matrix) to use with grofit
    clone_growth_data_time <- matrix(data = unique(clone_growth_data$Elapsed), 
                                     ncol = length(unique(clone_growth_data$Elapsed)), 
                                     nrow = nrow(clone_growth_data_rate), byrow = T)
    # run grofit
    clone_growth_data_grofit <- grofit(clone_growth_data_time, clone_growth_data_rate, control = grofit.control(interactive=F))
    
    # grofit summary
    clone_growth_data_grofitsum <- summary(clone_growth_data_grofit$gcFit) %>%
      select(1:9, 13)
    
    colnames(clone_growth_data_grofitsum)[c(1:3, 9, 10)] <- c('Plate', 'Well', 'Content', 'mu', 'stdmu') #mu is the maximum slope, and stdmu its std deviation
     # NOTE: sometimes the model can't be fitted and we get NA's in mu and stdmu, which are not plotted. 
     # Curves look ragged but fine, so we are losing data here due to lack of model fitting
    
  
    return(clone_growth_data_grofitsum)
}


# -----------------------------
# PLOTS

# clone growth curves
fun.growth_plotcurves <- funtion(cgd) {
  # cgd: clone_growth_data
  
clone_growth_data <- as.data.frame(cgd) %>%
  select(Plate, Well, Content, Elapsed, Confluence) %>%
  mutate(Content = as.factor(Content))

clone_growth_curve <- clone_growth_data %>% 
                      group_by(Content, Plate, Well) %>% #this is necessary to do the plotting. It gives something similiar to ggplot's aes(group=)
                      plot_ly(., 
                              x = ~Elapsed, 
                              y = ~Confluence,
                              color = ~Content,
                              type = 'scatter',
                              mode = 'lines',
                              # Hover text:
                              text = ~paste("Content:Plate:Well ", Content, Plate, Well)) %>%
                      layout(xaxis = list(title ='Time (h)'), yaxis = list(title = 'Confluence (%)'))

return(clone_growth_curve)
}

# clone growth rates
fun.growth_plotrates <- function(cgdg) {
  #cgdg: clone_growth_data_grofitsum

clone_growth_rate <- subset(clone_growth_data_grofitsum, !is.na(mu)) %>%
                     group_by(Content, Plate, Well) %>%
                     plot_ly(.,
                             x = ~Content,
                             y = ~mu,
                             type = 'scatter',
                             mode = 'markers',
                             # Hover text:
                             text = ~paste("Content:Plate:Well ", Content, Plate, Well)
                             ) %>%
                     layout(xaxis = list(title ='Clone'), yaxis = list(title = 'Maximum growth slope'))

return(clone_growth_rate)
}