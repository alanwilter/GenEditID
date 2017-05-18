# CRAN libraries
library(shiny)    # app
library(ggplot2)  # plotting
library(plotly)   # interactive plots over ggplots
library(reshape2)
library(dplyr)
library(RPostgreSQL)
library(DT)
library(grofit)   # calculation of growth slopes

# --------------------------------------------------------------------------------
# SCRIPTS. Loading of functions used for data processing and plotting

scriptsDir <- file.path(dirname(getwd()), "scripts")

source(file.path(scriptsDir, 'growth_plot.R'))
  # fun.growth_readDB
  # fun.growth_calcslope
  # fun.growth_plotcurves
  # fun.growth_plotrates

source(file.path(scriptsDir, 'abundance_plot.R'))
  # fun.protein_readDB
  # fun.protein_calcratio
  # fun.protein_plot

source(file.path(scriptsDir, 'NGS_plot.R'))
  # fun.NGS_readDB
  # fun.NGS_indelranges
  # fun.NGS_exploratory
  # fun.NGS_plotindelranges
  # fun.NGS_plotmutations
  # fun.NGS_plotvariants
  # fun.NGS_plotzygosity
  # fun.NGS_plotdistance

# source('layout_plot.R')

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
        
        # --------------------------------------------------------------------------------
        # PLOTS
        
        # plot clone growth curves
        clone_growth_data <- fun.growth_readDB(db)
        #clone_growth_data$cLayout <- as.factor(sapply(clone_growth_data$Content, function(a) strsplit(a, " ")[[1]][3]))
        clone_growth_curve <- fun.growth_plotcurves(clone_growth_data)
        
        # plot clone growth rates
        clone_growth_data_grofitsum <- fun.growth_calcslope(clone_growth_data)
        clone_growth_rate <- fun.growth_plotrates
        
        # plot ratio against cell line
        protein_abundance_data <- fun.protein_readDB(db) %>% fun.protein_calcratio()
        protein_abundance_plot <- fun.protein_plot(protein_abundance_data)
          
        #--------
        # NGS plots
        
        # indel positional ranges
        NGS.data <- fun.NGS_readDB(db)
        NGS.indel_ranges <- fun.NGS_indelranges(NGS.data)
        NGS.positional_plot <- fun.NGS_plotindelranges(NGS.indel_ranges)
        
        NGS.exploratory <- fun.NGS_exploratory(db)
        
        # exploratory plots: Type of mutation.
        #NGS.mutations_plot <- fun.NGS_plotmutations(---)
        
        # exploratory plots: Type of variant (frameshift, in-frame deletion, insertion, etc).
        #NGS.variants_plot <- fun.NGS_plotvariants(---)
        
        # exploratory plots: zygosity (homo, muthet, mutwt). 
        # Note: a muthet has two mutated alleles, each with a different mutation. A mutwt has one mutant allele and one wt allele
        #NGS.zygosity_plot <- fun.NGS_plotzygosity(---)
        
        # exploratory plots: distance of indel to cut site
        #NGS.distance_plot <- fun.NGS_plotdistance(---)
        
        
        # plot combined data (growth + NGS + protein)
        
        # plot combined data (growth + protein)
        
        # plot combined data (protein + NGS)
        
        # plot results layout plate
        
        
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
