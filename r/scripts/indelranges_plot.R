# DESCRIPTION:
  # R script to process NGS indel and SNV data.
  # Functions to read from database, process the data and and create a tile plot with indel positions and allele frequencies per well.

# Currently it does not read from the database. Also we need to modify the hard-added column 'Type', so it is managed from the database only.
# Things that should come or go from/to the database are marked as %DB

library(dplyr)
library(plotly)


fun.NGS_readDB <- function() {
  
  # ----------update this to read from DB
  
  # Hard-added a column called 'Type' with the variant caller used (Haplotype.caller or Vardict). This should be done differently when loading from
  # database, since data from indels and snvs is unified there
  NGSdata <- read.delim('r/shinyapp/data/SLX-13774variants.csv')
  NGSdata <- NGSdata[,c(1:4, 9, 10, 16,17, 20, 31)]
  colnames(NGSdata) <- c('Type', 'Sample', 'Barcode', 'Variant', 'Chromosome', 'Position', 'Allele.fraction', 'Depth', 'Amplicon', 'indel.length')
  # ----------
  
  return(data)
}

#-----------------------------
# CALCULATIONS

fun.NGS_calc <- function(NGSdata) {
  # NGSdata: data read from the database

    # ---------Ruben update this once the database table is clear  
    # take results from Haplotype.caller only
    data <-  droplevels(NGSdata[grep('H.', NGSdata$Type), ]) # split is done by levels, so make sure that unused ones are dropped
    # ---------
    
    databy <- split(data, f = data$Sample)
    
    # create indel position ranges (position, indellenght), shifting positions for deletions
     # (e.g. an indel in position 8 and width = -2 is converted into indel in postion 6 and width = 2) for IRanges to work
    databy.haploranges <- lapply(databy, function(b) {
      dat <- apply(b[c(6,10)], 1, FUN = function(a) {    #note that `apply` works on matrices, so it must receive just numeric data, otherwise it will convert 'b' into character.
        if(is.na(a[2])) res <- c(a[1], a[1], 0) else
          if(a[2] < 0) res <- c(a[1], a[1] + a[2], abs(a[2])) else
            if(a[2] > 0) res <- c(a[1], a[1], a[2])
            names(res) <- c('Original.position', 'Position.IRanges', 'indel.IRanges.width') # the 'original.position column is required for the 'pseudomerge' operation in line 64
            return(res)}
          ) %>% t  #to transpose the matrix
    })
    
    # identify each allele individually in each sample (*)
    databy.haploranges <- mapply(function(a,b) {a <- as.data.frame(a)
                                    a$indelID <- letters[1:nrow(a)] #(*)
                                    a$Sample <- b
                                    rownames(a) <- NULL
                                    return(a)},
                                 databy.haploranges, names(databy.haploranges),
                                 SIMPLIFY = F
                                 )
    
    databy.haploranges.nalleles <- sapply(databy.haploranges, function(a) nrow(a))  # number of alleles per sample
    databy.haploranges.ranges <-  do.call(rbind, databy.haploranges) %>%            # indel ranges
      mutate(start = Position.IRanges, end = Position.IRanges + indel.IRanges.width)
    
    # merge position ranges with original indel data
     # note that sample names are duplicated (because there are several indels per sample), 
     # so `merge` or `inner_join` will fail and generate incorrect (bigger) dataframes unless more index columns are used to make rows unique.
     # Instead of that, merge by ordering and binding the data.
    databy.haploranges.ranges <- databy.haploranges.ranges[order(databy.haploranges.ranges$Sample, databy.haploranges.ranges$Original.position), ]
    data <- data[order(data$Sample, data$Position),]
    databy.haploranges.rangesfull <- data.frame(data, databy.haploranges.ranges)
    
    # delete columns 'Sample.1' and 'Original.position'
    databy.haploranges.rangesfull <- databy.haploranges.rangesfull[, -grep('Sample.1|Original.position', colnames(databy.haploranges.rangesfull)),]

    return(databy.haploranges.rangesfull)
}

# -----------------------------
# PLOTS

fun.NGS_plotindelranges <- function(a) {
  # a: databy.haploranges.rangesfull

# ------ Ruben to update ggplot with plot_ly  
  g.plot <- plot_ly(subset(a, indel.IRanges.width > 0),
                    shapes = list(
                      list(
                           type = "rect",
                           fillcolor = "blue", line = list(color = "blue"), opacity = 0.3,
                           x0 = "1980-01-01", x1 = "1985-01-01", xref = "x",
                           y0 = 4, y1 = 12.5, yref = "y"),
                      )
                    )
  
  
  
g.plot <- ggplot(subset(a, indel.IRanges.width > 0), aes(xmin = start, xmax = end, ymin = 0, ymax = Allele.fraction, color = indelID)) +
  geom_point(data = subset(a, indel.IRanges.width == 0 & Chromosome == 'chr10'), aes(x = Position, y = Allele.fraction), color = 'black') +
  geom_rect(mapping = aes(fill = Variant)) +
  facet_wrap(~Sample) +
  geom_vline(xintercept=89653802, color = 'steelblue', linetype = "dotted") +
  theme_classic()

  return(g.plot)
}