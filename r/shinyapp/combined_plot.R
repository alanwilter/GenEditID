# DESCRIPTION:
# R script to combine growth, protein and NGS data in pairs or in a three-variable plot
# Functions to read from database, process the data and:

######### not finished !!!------------

library(dplyr)
library(plotly)


# Data is loaded through glowth_plot.R, abundance_plot.R and NGS_plot.R

#-----------------------------
# CALCULATIONS
# merge ICW and incu data in a single table for plotting:
clone_growth_protein <- merge(clone_growth_data_grofitsum, protein_abundance_data, all = T, by = c('Plate', 'Well', 'Content')) %>%
  mutate('minusSD' = mu - stdmu, 'plusSD' = mu + stdmu) %>%
  subset(!grepl('normalisation|background', Content))

# -----------------------------
# PLOTS

# # combined plot (protein + growth slopes + NGS)
# ICWincuNGS_plotly <- plot_ly(ICWincuNGS, type = 'scatter', mode = 'markers',
#                              # Hover text:
#                              text = ~paste("Plate:Well ", Plate, Well)) %>%
#   add_trace(x = ~ratio800to700, y = ~mu, symbol = ~Content, color = I("black")) %>% #error bars not showing in the right place... error_y = ~list(type = "data", array = stdmu)
#   add_trace(subset(ICWincuNGS, !is.na(Zygosity)), x = ~ratio800to700, y= ~mu, color = ~Zygosity,
#             text = ~paste("Plate:Well:Zygosity ", Plate, Well, Zygosity)) %>%
#   layout(xaxis = list(title ='protein abundance'), yaxis = list(title = 'maximum growth slope'))
# 
