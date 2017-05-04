library(dplyr)
library(grofit)
library(ggplot2)

# connect to database
#db <- src_sqlite("crispr.sqlite")
db <- src_postgres(user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")

# connect to tables
cell_line <- tbl(db, "cell_line") %>% rename(cell_line_id=id, cell_line_name=name, cell_line_description=description)
project <- tbl(db, "project") %>% rename(project_id=id, project_geid=geid, project_description=description)
experiment_layout <- tbl(db, "experiment_layout") %>% rename(experiment_layout_id=id, experiment_layout_geid=geid)
clone <- tbl(db, "clone") %>% rename(clone_id=id, clone_name=name)
well <- tbl(db, "well") %>% rename(well_id=id)
well_content <- tbl(db, "well_content") %>% rename(well_content_id=id)
plate <- tbl(db, "plate") %>% rename(plate_id=id)
growth <- tbl(db, "growth") %>% rename(growth_id=id)
abundance <- tbl(db, "abundance") %>% rename(abundance_id=id)
target <- tbl(db, "target") %>% rename(target_id=id, target_name=name)
guide <- tbl(db, "guide") %>% rename(guide_id=id, guide_name=name)
guide_well_content_association <- tbl(db, "guide_well_content_association")

create_classifier <- function(cellline, clone, guide, content)
{
    thing <- paste(cellline, clone, guide, content, sep=" ")
    thing <- gsub(" NA", "", thing)
    thing <- gsub(" sample", "", thing)
    thing
}

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
    mutate(position=paste0(row, column)) %>%
    mutate(classifier=create_classifier(cell_line_name, clone_name, guide_name, content_type)) %>%
    select(classifier, position, intensity_channel_700, intensity_channel_800) %>%
    rename(Content=classifier, Channel700=intensity_channel_700, Channel800=intensity_channel_800)

print(protein_abundance_data)

# calculate relative protein abundance, ratio 800 to 700
protein_abundance_data['ratio800to700'] = protein_abundance_data['Channel800']/protein_abundance_data['Channel700']
protein_abundance_data

# plot ratio against cell line
protein_abundance_plot <- ggplot(protein_abundance_data, aes(x = Content, y = ratio800to700)) +
    geom_violin() +
    geom_point(aes(colour=Content)) +
    theme_bw() +
    ggtitle('InCellWestern') +
    xlab('Cell line') +
    ylab('Relative protein abundance (ratio 800 to 700 nm)') +
    theme(axis.text.x = element_text(angle=15, vjust=0.5, size=10))

print(protein_abundance_plot)

# read clone growth data
clone_growth_data <- left_join(well, growth, by="well_id") %>%
    inner_join(., well_content, by="well_content_id") %>%
    inner_join(., clone, by="clone_id") %>%
    inner_join(., plate, by="plate_id") %>%
    inner_join(., guide_well_content_association, by="well_content_id") %>%
    inner_join(., guide, by="guide_id") %>%
    inner_join(., target, by="target_id") %>%
    collect %>%
    mutate(position=paste0(row, column)) %>%
    select(plate_id, target_name, guide_name, position, hours, confluence_percentage) %>%
    rename(PlateID=plate_id, Target=target_name, Guide=guide_name, Well=position, Elapsed=hours, Confluence=confluence_percentage) %>%
    filter(!is.na(Confluence))

clone_growth_data

# plot clone growth curves
clone_growth_curve <- ggplot(clone_growth_data, aes(x = Elapsed, y = Confluence, group = as.factor(Target):as.factor(PlateID):as.factor(Well))) +
    geom_line(aes(colour = Target)) +
    theme_bw() +
    ggtitle('Clone growth curves') +
    xlab('Elapse time') +
    ylab('Confluence (%)')

print(clone_growth_curve)

# TODO create one plot per well with plotly

# calculate slope of growth curves with grofit
# Moving average: https://en.wikipedia.org/wiki/Moving_average
# Kernel smoother: https://en.wikipedia.org/wiki/Kernel_smoother
# TODO

# plot clone growth rate
# TODO
