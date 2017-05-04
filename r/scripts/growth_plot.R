library(dplyr)
library(grofit)
library(ggplot2)

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
