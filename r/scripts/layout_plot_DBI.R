library(dplyr)
library(ggplot2)
library(svglite)
library(RColorBrewer)
library(RPostgreSQL)


## Plot function, stolen from "PlateLayout".

testPlot <- function(plotDat, plateName, palette, saveTo = NULL)
{
    conditionToPlot <- "Content"

    ContentFactors = as.factor(plotDat$content)
    plotDat$Content <- factor(ContentFactors, levels=levels(ContentFactors))
    
    boxColCodes <- colorRampPalette(palette)(length(levels(ContentFactors)))
    
    RowLabels <- as.factor(plotDat$row)
    plotDat$row <- factor(RowLabels, levels=rev(levels(RowLabels)))
    
    NumberOfColumns <- max(plotDat$column)

    passes <- plotDat$ge_score >= 3
    boxColCodes <- factor(ifelse(passes, "green", "red"))
    
    print(boxColCodes)
    
    plot <- ggplot(plotDat) +
        geom_tile(aes(x=column, y=row, fill=passes), colour="black")  +
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
    thing <- gsub(" empty", "", thing)
    thing
}


loadPlateData <- function(conn, layoutId)
{
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


## Main code.

driver <- dbDriver("PostgreSQL")

conn <- dbConnect(driver, user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")

layoutId <- 'GEP00001_03'

tryCatch(
    {
        plateDat <- loadPlateData(conn, layoutId)
        
        plot <- testPlot(plateDat, layoutId, brewer.pal(12, "Set3"), paste0(layoutId, "_layout.svg"))
        
        print(plot)
    },
    finally =
    {
        dbDisconnect(conn)
    }
)
