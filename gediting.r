library(dplyr)
library(grofit)
library(ggplot2)

# connect to database
db <- src_sqlite("crispr.sqlite")

# connect to tables
clone <- tbl(db, "clone") %>% 
  rename(clone_id=id)
well <- tbl(db, "well") %>% 
  rename(well_id=id)
growth <- tbl(db, "incucyte") 

# read protein abundance data
protein_abundance_data <- left_join(well, clone, by="clone_id") %>% 
  collect %>% 
  mutate(position=paste0(row, column)) %>% 
  select(plate_id, name, position, icw_700, icw_800) %>% 
  rename(PlateID=plate_id, Target=name, Well=position, Channel700=icw_700, Channel800=icw_800)

print(protein_abundance_data)

# calculate relative protein abundance, ratio 800 to 700
protein_abundance_data['ratio800to700'] = protein_abundance_data['Channel800']/protein_abundance_data['Channel700']
protein_abundance_data

# plot ratio against cell line
protein_abundance_plot <- ggplot(protein_abundance_data, aes(x = Target, y = ratio800to700)) +
  geom_violin() +
  geom_point(aes(colour=Target)) +
  theme_bw() + 
  ggtitle('InCellWestern') +
  xlab('Cell line') +
  ylab('Relative protein abundance (ratio 800 to 700 nm)') +
  theme(axis.text.x  = element_text(angle=15, vjust=0.5, size=10))

print(protein_abundance_plot)

# read clone growth data
clone_growth_data <- left_join(well, clone, by="clone_id") %>% 
  left_join(growth, well, by="well_id") %>% 
  collect %>% 
  mutate(position=paste0(row, column)) %>% 
  select(plate_id, name, position, hours, phased_object_confluence) %>% 
  rename(PlateID=plate_id, Target=name, Well=position, Elapsed=hours, Confluence=phased_object_confluence) %>% 
  filter(!is.na(Confluence))
  
clone_growth_data 

# plot clone growth curves
clone_growth_curve <-  ggplot(clone_growth_data, aes(x = Elapsed, y = Confluence, group = as.factor(Target):as.factor(PlateID):as.factor(Well))) +
  geom_line(aes(colour = Target)) +
  theme_bw() +
  ggtitle('Clone growth curves') +
  xlab('Elapse time') +
  ylab('Confluence (%)')

print(clone_growth_curve)

# TODO create one plot per well with plotly

# calculate slope of growth curves with grofit
# TODO

# plot clone growth rate
# TODO
