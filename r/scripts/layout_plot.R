library(dplyr)
library(ggplot2)
library(svglite)
library(RColorBrewer)

## Plot function, stolen from "PlateLayout".

testPlot <- function(plotDat, plateName, palette, saveTo = NULL)
{
    conditionToPlot <- "Content"

    ContentFactors = as.factor(plotDat$Content)
    plotDat$Content <- factor(ContentFactors, levels=levels(ContentFactors))
    
    boxColCodes <- colorRampPalette(palette)(length(levels(ContentFactors)))
    
    RowLabels <- as.factor(plotDat$Row)
    plotDat$Row <- factor(RowLabels, levels=rev(levels(RowLabels)))
    
    NumberOfColumns <- max(plotDat$Column)

    passes <- plotDat$Score >= 3
    boxColCodes <- factor(ifelse(passes, "green", "red"))
    
    print(boxColCodes)
    
    plot <- ggplot(plotDat) +
        geom_tile(aes(x=Column, y=Row, fill=passes), colour="black")  +
        scale_x_continuous(breaks=1:NumberOfColumns, position="top") +
        scale_fill_manual(values=levels(boxColCodes), na.translate=FALSE, labels=c("<3", ">=3")) +
        theme(axis.text.x=element_text(size=20), axis.text.y=element_text(size=20),
              axis.title.x=element_blank(), axis.title.y=element_blank()) +
        ggtitle(paste0(plateName, " Plate Layout")) +
        labs(fill="Score")
    
    if (!is.null(saveTo))
    {
        plotWidth <- 1.1 * NumberOfColumns + 1.5
        ggsave(saveTo, height = 5, width = plotWidth)
    }
    
    plot
}


create_classifier <- function(cellline, clone, guide, content)
{
    thing <- paste(cellline, clone, guide, content, sep=" ")
    thing <- gsub(" NA", "", thing)
    thing <- gsub(" sample", "", thing)
    thing
}


loadPlateData <- function(db, plateId)
{
    left_join(well, abundance, by="well_id") %>%
    left_join(., well_content, by="well_content_id") %>%
    left_join(., clone, by="clone_id") %>%
    left_join(., cell_line, by="cell_line_id") %>%
    left_join(., plate, by="plate_id") %>%
    left_join(., guide_well_content_association, by="well_content_id") %>%
    left_join(., guide, by="guide_id") %>%
    left_join(., target, by="target_id") %>%
    left_join(., sequencing_library_content, by="well_id") %>%
    left_join(., variant_result, by="sequencing_library_content_id") %>%
    collect %>%
    filter(plate_geid == plateId) %>%
    filter(!is.na(content_type)) %>%
    mutate(classifier=create_classifier(cell_line_name, clone_name, guide_name, content_type)) %>%
    select(classifier, row, column, ge_score) %>%
    rename(Content=classifier, Row=row, Column=column, Score=ge_score)
}


## Main code.

db <- src_postgres(user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")

cell_line <- tbl(db, "cell_line") %>% rename(cell_line_id=id, cell_line_name=name, cell_line_description=description)
project <- tbl(db, "project") %>% rename(project_id=id, project_geid=geid, project_description=description)
experiment_layout <- tbl(db, "experiment_layout") %>% rename(experiment_layout_id=id, experiment_layout_geid=geid)
clone <- tbl(db, "clone") %>% rename(clone_id=id, clone_name=name)
well <- tbl(db, "well") %>% rename(well_id=id)
well_content <- tbl(db, "well_content") %>% rename(well_content_id=id)
plate <- tbl(db, "plate") %>% rename(plate_id=id, plate_geid=geid)
growth <- tbl(db, "growth") %>% rename(growth_id=id)
abundance <- tbl(db, "abundance") %>% rename(abundance_id=id)
target <- tbl(db, "target") %>% rename(target_id=id, target_name=name)
guide <- tbl(db, "guide") %>% rename(guide_id=id, guide_name=name)
guide_well_content_association <- tbl(db, "guide_well_content_association")
sequencing_library_content <- tbl(db, "sequencing_library_content") %>% rename(sequencing_library_content_id=id)
variant_result <- tbl(db, "variant_result") %>% rename(variant_result_id=id)


## TODO: Pass plate identifier as argument to script.

plateId <- 'GEP00001_03'

plateDat <- loadPlateData(db, plateId)

plot <- testPlot(plateDat, plateId, brewer.pal(12, "Set3"), paste0(plateId, "_layout.svg"))

print(plot)

RPostgreSQL::dbDisconnect(db$con)
