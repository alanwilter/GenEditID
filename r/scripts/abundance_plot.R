library(dplyr)
library(grofit)
library(ggplot2)

# connect to database
db <- src_postgres(user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")

# connect to tables
cell_line <- tbl(db, "cell_line") %>% rename(cell_line_id=id, cell_line_name=name, cell_line_description=description)
clone <- tbl(db, "clone") %>% rename(clone_id=id, clone_name=name)
well <- tbl(db, "well") %>% rename(well_id=id)
well_content <- tbl(db, "well_content") %>% rename(well_content_id=id)
plate <- tbl(db, "plate") %>% rename(plate_id=id)
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
