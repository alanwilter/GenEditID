# DESCRIPTION:
  # R script to process protein abundance data from the Incell Western.
  # Functions to read from database, calculate protein relative abundance (ratio 800 to 700 nm) and produce a plots of protein abundance.

library(dplyr)
library(plotly)

fun.protein_readDB <- function(db) {
  # connect to tables
  well <- tbl(db, "well") %>% rename(well_id=id)
  abundance <- tbl(db, "abundance") %>% rename(abundance_id=id)
  well_content <- tbl(db, "well_content") %>% rename(well_content_id=id)
  clone <- tbl(db, "clone") %>% rename(clone_id=id, clone_name=name)
  cell_line <- tbl(db, "cell_line") %>% rename(cell_line_id=id, cell_line_name=name, cell_line_description=description)
  plate <- tbl(db, "plate") %>% rename(plate_id=id)
  guide_well_content_association <- tbl(db, "guide_well_content_association")
  guide <- tbl(db, "guide") %>% rename(guide_id=id, guide_name=name)
  target <- tbl(db, "target") %>% rename(target_id=id, target_name=name)
  
  # read protein abundance data
  protein_abundance_data <- left_join(well, abundance, by="well_id") %>%
    left_join(., well_content, by="well_content_id") %>%
    left_join(., clone, by="clone_id") %>%
    left_join(., cell_line, by="cell_line_id") %>%
    left_join(., plate, by="plate_id") %>%
    left_join(., guide_well_content_association, by="well_content_id") %>%
    left_join(., guide, by="guide_id") %>%
    left_join(., target, by="target_id") %>%
    collect %>%
    filter(!is.na(content_type)) %>%
    mutate(Well = paste0(row, column)) %>%
    mutate(classifier=create_classifier(cell_line_name, clone_name, guide_name, content_type)) %>%
    select(classifier, plate_id, Well, intensity_channel_700, intensity_channel_800) %>%
    rename(Content=classifier, Plate = plate_id, Channel700=intensity_channel_700, Channel800=intensity_channel_800) %>%
    filter(!grepl('background', Content)) %>%            #the background control (no antibodies is not meaningful for proteini amount)
    mutate('ratio800to700' = Channel800/Channel700) %>%  #calculate relative protein abundance, ratio 800 to 700
    mutate(Plate = as.factor(Plate), Well = as.factor(Well))  #change Plate and position to factors, so grouping in ggplot can work
  
  return(protein_abundance_data)
}

#-----------------------------
# CALCULATIONS

# calculate relative protein abundance, ratio 800 to 700
fun.protein_calcratio <- function(pad) {
  # pad: protein_abundance_data
  
    pad['ratio800to700'] <- pad['Channel800']/pad['Channel700']
  return(pad)
}

# -----------------------------
# PLOTS

# plot ratio against cell line
fun.protein_plot <- function(pad) {
  # pad: protein_abundance_data
  
  protein_abundance_plot <- pad %>%
                          group_by(Content, Plate, Well) %>%
                          plot_ly(.,
                                  x = ~Content,
                                  y = ~ratio800to700,
                                  type = 'scatter',
                                  mode = 'markers',
                                  color = ~Content,
                                  # Hover text:
                                  text = ~paste("Content:Plate:Well ", Content, Plate, Well)
                                  ) %>%
                          layout(xaxis = list(title ='Cell line', showticklabels = FALSE), yaxis = list(title = 'Relative protein abundance'))

return(protein_abundance_plot)
}