# CRAN libraries
library(shiny)    # app
library(ggplot2)  # plotting
library(plotly)   # interactive plots over ggplots
library(grid)     # used in NGSplots.R, neccesary for function grid::grid.draw
library(gridExtra)# for arrangement of ggplots in a grid
library(reshape2)
library(dplyr)
library(RPostgreSQL)
library(DT)
library(grofit)   # calculation of growth slopes
library(DBI)      # database interface definition for communication between R and relational database management systems
library(RColorBrewer) # to create color palettes for plotting

# --------------------------------------------------------------------------------
# SCRIPTS. Loading of functions used for data processing and plotting

#scriptsDir <- file.path(dirname(getwd()), "scripts") # if in shiny, all paths are relative to the shiny root directory

#source(file.path(scriptsDir, 'growth_plot.R'))
source('growth_plot.R')
  # fun.growth_readDB
  # fun.growth_calcslope
  # fun.growth_plotcurves
  # fun.growth_plotrates

#source(file.path(scriptsDir, 'abundance_plot.R'))
source('abundance_plot.R')
  # fun.protein_readDB
  # fun.protein_calcratio
  # fun.protein_plot

#source(file.path(scriptsDir, 'NGS_plot.R'))
source('NGS_plot.R')
  # fun.NGS_readDB
  # fun.NGS_indelranges
  # fun.NGS_exploratory
  # fun.NGS_plotindelranges
  # fun.NGS_plotmutations
  # fun.NGS_plotvariants
  # fun.NGS_plotzygosity
  # fun.NGS_plotdistance
  # fun.NGS_plotalleles

source('layout_plot.R')
  # fun.layoutplot

# source('combined_plot.R')

# --------------------------------------------------------------------------------
# FUNCTIONS

create_classifier <- function(cellline, clone, guide, content)
{
  thing <- paste(cellline, clone, guide, content, sep=" ")
  thing <- gsub(" NA", "", thing)
  thing <- gsub(" sample", "", thing)
  thing <- gsub(" empty", "", thing)
  thing
}

# Database connection.
# Note that this needs to be safely released, hence tryCatch.

db <- src_postgres(user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")

tryCatch(
    {
        # ------------------
        # DATA processing
        
        # NGS data
          # merge data from funNGS.indelranges and fun.NGS_exploratory. Create 'ids' factor for merging in NGSdata
        NGSdata <- fun.NGS_readDB(db)
        NGSdata.full <- mutate(NGSdata, 'ids' = do.call(paste, c(list(sample, type, guide), list(sep = '.')))) %>%
          inner_join(., fun.NGS_indelranges(NGSdata), by = c('sample', 'position')) %>%
          inner_join(., fun.NGS_exploratory(NGSdata), by = c('ids', 'gene'))
        
        
        # ------------------
        # PLOTS
        
        # plot clone growth curves
        clone_growth_data <- fun.growth_readDB(db)
        #clone_growth_data$cLayout <- as.factor(sapply(clone_growth_data$Content, function(a) strsplit(a, " ")[[1]][3]))
        clone_growth_curve <- fun.growth_plotcurves(clone_growth_data)
        
        # plot clone growth rates
  #      clone_growth_data_grofitsum <- fun.growth_calcslope(clone_growth_data)
   #     clone_growth_rate <- fun.growth_plotrates
        
        # plot ratio against cell line
        protein_abundance_data <- fun.protein_readDB(db) %>% fun.protein_calcratio()
        protein_abundance_plot <- fun.protein_plot(protein_abundance_data)
          
        #--------
        # NGS plots
        
        ## exploratory plots
        NGS.plotindel_ranges <- fun.NGS_plotindelranges(NGSdata.full) #note that sql query for this needs fixing (see planning.md)
        NGS.plotmutations    <- fun.NGS_plotmutations(NGSdata.full)
        NGS.plotvariants     <- fun.NGS_plotvariants(NGSdata.full)
        NGS.plotzygosity     <- fun.NGS_plotzygosity(NGSdata.full)
        NGS.plotdistance     <- fun.NGS_plotdistance(NGSdata.full)
        NGS.plotalleles      <- fun.NGS_plotalleles(NGSdata.full)
        
        
        # plot combined data (growth + NGS + protein)
        
        # plot combined data (growth + protein)
        
        # plot combined data (protein + NGS)
        
        # plot results layout plate
        layoutId <- 'GEP00001_03' # this should extend to all plates in db
        plate.plotscores <- fun.layoutplot(conn, layoutId)
        
        # --------------------------------------------------------------------------------
        # Notes
        # TODO create one plot per well with plotly
        # TODO calculate slope of growth curves with grofit
        # TODO Moving average: https://en.wikipedia.org/wiki/Moving_average
        # TODO Kernel smoother: https://en.wikipedia.org/wiki/Kernel_smoother
    },
    finally =
    {
        RPostgreSQL::dbDisconnect(db$con)
    }
)
