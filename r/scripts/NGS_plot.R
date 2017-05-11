# DESCRIPTION:
  # R script to process NGS indel and SNV data.
  # Functions to read from database, process the data and:
  # 1. create a tile plot with indel positions and allele frequencies per well.
  # 2. calculate types of mutations (homo, muthet, mutwt).
  # 3. make exploratory graphs

library(dplyr)
library(plotly)

# Currently it does not read from the database. Also we need to modify the hard-added column 'Type', so it is managed from the database only.
# Things that should come or go from/to the database are marked as %DB

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

# calculate indel positional ranges for fun.NGS_plotindelranges
fun.NGS_indelranges <- function(NGSdata) {
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

# calculate data for exploratory plots
fun.NGS_exploratory <- function(NGSdata) {
  # NGSdata: data read from the database
  
  
  # Should be loaded from DB
  # Things annotated as %DB% and /%DB% require updating after NGS data is in the database
  # %DB%
  # Temporary workaround to get the NGS data, till we have it in the database.
  # To create the RDS files, look into r/scripts/ICWincuNGSdata/incu_incellNGS_GEP00001.R script to recreate them
  data <- readRDS("data/GEP00001_data.RDS")
  data <- data[data$Allele.fraction > 0.15, c('Plate', 'Well', 'guide', 'Alleles', 'Indel.length',
                                              'Allele.fraction', 'Chromosome', 'Position',
                                              'Variant.type.consequence', 'Symbol..Gene.ID.')]
  data$Content <- paste0('MCF7 clone3 ', data$guide)
  data <- data[order(data$Plate, data$Well), ]
  
  NGSdata <- readRDS('data/GEP00001_dataNGS.RDS') %>%
    mutate('Content' = as.factor(paste0('MCF7 clone3 ', Layout)),
           'Plate' = as.factor(Plate)) %>%
    mutate('Plate' = gsub('plate', '', Plate)) %>%
    droplevels %>%
    arrange(desc(Plate, Well))
  
  NGSdata.cells <- NGSdata[grepl('-C', NGSdata$Sample),] # extracted cells only
  NGSdata.gDNA <- NGSdata[grepl('-G', NGSdata$Sample),]  # extracted gDNA only
  
  fun.by <- function(x) {NGS <- by(x[,], INDICES = list(x$Plate, x$Well, x$Type), FUN = as.data.frame)
  NGS <- NGS[!do.call(c, lapply(NGS, is.null))]
  NGS <- lapply(NGS, function(a) data.frame(
    'Plate' = a$Plate,
    'Well' = gsub('P.', '', a$Well),
    'Content' = paste0('MCF7 clone3 ', a$Layout),
    'Type' = a$Type,  #SNV or INDEL
    'NROW' = nrow(a), 
    'homo' = a$AlleleFraction >0.85,
    'muthet' = a$AlleleFraction > 0.35 & a$AlleleFraction <0.85 & nrow(a) > 1,
    'mutwt' = a$AlleleFraction > 0.35 & a$AlleleFraction <0.85 & nrow(a) == 1,
    'has.offtargets' = length(unique(a$Gene)) > 1)
  )
  NGS2 <- lapply(NGS, function(a) {  #consolidate rows of the dataframe
    a$homo <- all(a$homo)
    a$muthet <- as.logical(sum(a$muthet))   #if it has three rows, two T and one F (e.g. it's a muthet with an off target), it will show as muthet
    a$mutwt <- as.logical(sum(a$mutwt))
    a$has.offtargets <- all(a$has.offtargets)
    return(a)
  })
  NGS <- lapply(NGS, function(a) a[1,]) #select first row only
  }
  
  fun.by2 <- function(x) {NGS <- by(x[,], INDICES = list(x$Plate, x$Well), FUN = as.data.frame)
  NGS <- NGS[!do.call(c, lapply(NGS, is.null))]}
  
  NGSdataby.cells <- lapply(fun.by(NGSdata.cells), function(a) data.frame(a, 'DNAsource' = 'cells'))
  NGSdataby.gDNA <- lapply(fun.by(NGSdata.gDNA), function(a) data.frame(a, 'DNAsource' = 'gDNA'))
  NGSdataby <- rbind(do.call(rbind, NGSdataby.cells),
                     do.call(rbind, NGSdataby.gDNA))
  NGSdataby$Zygosity <- apply(NGSdataby[,c('homo', 'muthet', 'mutwt')], 1, function(a) c('homo', 'muthet', 'mutwt')[a]) %>%
    lapply(function(a) {if(length(a) == 0) a <- NA; return(a)}) %>%
    do.call(rbind, .) %>%
    as.factor
  
  # to check odd results
  # m <- NGSdataby[duplicated(paste0(NGSdataby$Plate, NGSdataby$Well)), c('Plate', 'Well')]
  # m <- paste0(m$Plate, m$Well)
  # NGSdataby$m2 <- as.character(paste0(NGSdataby$Plate, NGSdataby$Well))
  # subset(NGSdataby, m2 %in% m) %>% arrange(Plate, Well) %>% group_by(Plate, Well)
  # subset(NGSdata.cells, Well == 'P3D7')
  # subset(NGSdataby, m2 == '3D7') %>% arrange(Plate, Well) %>% group_by(Plate, Well)
  
  # NGSdataby.cells[do.call(rbind, lapply(NGSdataby.cells, function(a) nrow(a) > 2))]
  # fun.by2(NGSdata.cells)[do.call(rbind, lapply(NGSdataby.cells, function(a) nrow(a) > 2))]
  
  # ICWincuNGS <- merge(clone_growth_protein, NGSdataby, by = c('Plate', 'Well', 'Content'), all = T)
  # ICWincuNGS$plusSD <- ICWincuNGS$mu + ICWincuNGS$stdmu
  # ICWincuNGS$minusSD <- ICWincuNGS$mu - ICWincuNGS$stdmu
  # # remove background and normalisation controls
  # ICWincuNGS <- subset(ICWincuNGS, !grepl('normalisation|background', Content)) %>% droplevels
  # #NGSdata <- subset(NGSdata, !grepl('normalisation|background', Content))
  # #/%DB%
  
}


# -----------------------------
# PLOTS

# indel positional ranges
fun.NGS_plotindelranges <- function(a) {
  # a: databy.haploranges.rangesfull

# ------ Ruben to potentially update ggplot with plot_ly  
          # g.plot <- plot_ly(subset(a, indel.IRanges.width > 0),
          #                   shapes = list(
          #                     list(
          #                          type = "rect",
          #                          fillcolor = "blue", line = list(color = "blue"), opacity = 0.3,
          #                          x0 = "1980-01-01", x1 = "1985-01-01", xref = "x",
          #                          y0 = 4, y1 = 12.5, yref = "y")
          #                     )
          #                   )
          # 
  
  
g.plot <- ggplot(subset(a, indel.IRanges.width > 0), aes(xmin = start, xmax = end, ymin = 0, ymax = Allele.fraction, color = indelID)) +
  geom_point(data = subset(a, indel.IRanges.width == 0 & Chromosome == 'chr10'), aes(x = Position, y = Allele.fraction), color = 'black') +
  geom_rect(mapping = aes(fill = Variant)) +
  facet_wrap(~Sample) +
  geom_vline(xintercept=89653802, color = 'steelblue', linetype = "dotted") +
  theme_classic()

  return(g.plot)
}

# exploratory plots: Type of mutation.
fun.NGS_plotmutations <- function(a) {
  # a: NGSdata.cells
  
  NGS.mutations <- plot_ly(a, type = 'histogram', x = ~Type) %>%
                   layout(xaxis = list(title = 'Type of mutation'), yaxis = list(title = '% of all mutations'))
}

# exploratory plots: Type of variant (frameshift, in-frame deletion, insertion, etc).
fun.NGS_plotvariants <- function(a) {
  # a: NGSdata.cells
  
  NGS.variants <-  plot_ly(mutate(a, Variant = c('frameshift', 'synonymous', 'inframe.del', 'intron.noncoding', 'stop.frameshift')[
  match(a$Variant, c('frameshift', 'synonymous','inframe_deletion', 'intron,non_coding_transcript', 'stop_gained,frameshift'))]),
  type = 'histogram', x = ~Variant) %>%
  layout(xaxis = list(title = 'Type of variant'),
         yaxis = list(title = '% of all variants'),
         margin = list(l = 150))
}

# exploratory plots: zygosity (homo, muthet, mutwt). 
  # Note: a muthet has two mutated alleles, each with a different mutation. A mutwt has one mutant allele and one wt allele

fun.NGS_plotzygosity <- function(a) {
  #a: NGSdataby
  
  NGS.zygosity <- plot_ly(mutate(a, has.offtargets = c('in-target', 'off-target')[match(a$has.offtargets, c(F, T))]),
                        type = 'histogram', x = ~Zygosity, color = ~has.offtargets) %>%
                        layout(xaxis = list(title = 'Zygosity'), yaxis = list(title = '%'))
}

# exploratory plots: distance of indel to cut site
fun.NGS_plotdistance <- function(a) {
  #a: data
  
NGS.distancetocutsite <- plot_ly(a, type = 'histogram', x = ~Indel.length) %>%
                         layout(xaxis = list(title = 'Indel length (nt)'), yaxis = list(title = '%'))
}
