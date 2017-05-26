# DESCRIPTION:
  # R script to process NGS indel and SNV data.
  # Functions to read from database, process the data and:
  # 1. create a tile plot with indel positions and allele frequencies per well.
  # 2. calculate types of mutations (homo, muthet, mutwt).
  # 3. make exploratory graphs

library(dplyr)
library(plotly)
library(DBI)
library(grid)
library(gridExtra)

# read NGS data from the database
fun.NGS_readDB <- function(db) {

  sql <- 
    'select
             p.geid as Plate,
             concat(w.row, w.column) as Well,
             g.name as Guide,
             slc.sequencing_sample_name as Sample,
             slc.sequencing_barcode as Barcode,
             vr.consequence as Consequence,
             vr.gene_id as SymbolGene,
             vr.cdna_effect as CDNAEffect,
             vr.protein_effect as ProteinEffect,
             vr.codons as Codons,
             vr.chromosome as Chromosome,
             vr.position as Position,
             vr.sequence_ref as Ref,
             vr.sequence_alt as Alt,
             vr.allele_fraction as AlleleFraction,
             vr.depth as Depth,
             vr.quality as Quality,
             vr.amplicon as Amplicon,
             vr.gene as Gene,
             vr.exon as Exon,
             vr.offset_from_primer_end as PrimerEndOffset,
             vr.indel_length as IndelLength,
             vr.sequence_alleles as Alleles,
             vr.variant_type as Type
         from
             well w
             left join sequencing_library_content slc on slc.well_id = w.id
             inner join variant_result vr on vr.sequencing_library_content_id = slc.id
             inner join experiment_layout el on w.experiment_layout_id = el.id
             left join plate p on p.experiment_layout_id = el.id
             inner join well_content wc on w.well_content_id = wc.id
             left join guide_well_content_association gwca on gwca.well_content_id = wc.id
             inner join guide g on gwca.guide_id = g.id
         where
             vr.allele_fraction > 0.1'
         
         # rows with allelefraction < 0.1 are likely PCR or sequencing errors
  
    NGSdata <- dbGetQuery(db$con, sql) 
  
    # TO FIX!
    # Plate should be filtered to yield only "NGS plates"
    NGSdata <- filter(NGSdata, grepl('incu', plate)) # temporary workaround until we fix the plate_NGS problems
    data.guide_location <- data.frame('guide' = unique(NGSdata$guide), 'guide_location' = c(40497619, 40497633, 40498696, 40497587), stringsAsFactors = F)
    NGSdata <- inner_join(NGSdata, data.guide_location, by = 'guide')

    return(NGSdata)
}


#-----------------------------
# CALCULATIONS

# calculate indel positional ranges for fun.NGS_plotindelranges
  # in the original data, deletions are in the format (position, -lenght), but that's not ideal for plotting and it fails if we convert
  # the data to IRanges format. This function converts that format into (position - length, length)

fun.NGS_indelranges <- function(NGSdata) {
  # NGSdata: data read from the database, produced by fun.NGS_readDB(db)

    NGSdata.bysample <- split(NGSdata, f = NGSdata$sample)
    
    # convert indel position ranges (position, indellength) to (position - deletionlength, deletion length) only when length is negative
    # (e.g. an indel in position 8 and width = -2 is converted into indel in postion 6 and width = 2) for IRanges to work

    NGSdata.indelranges <- lapply(NGSdata.bysample, function(b) {
      dat <- apply(b[c('position','indellength')], 1, FUN = function(a) {    #note that `apply` works on matrices, so it must receive just numeric data, otherwise it will convert 'b' into character.
        if(is.na(a[2])) res <- c(a[1], a[1], 0) else
          if(a[2] < 0) res <- c(a[1], a[1] + a[2], abs(a[2])) else
            if(a[2] > 0) res <- c(a[1], a[1], a[2])
            names(res) <- c('position', 'Position.IRanges', 'indel.IRanges.width') # the 'original.position column is required for the 'inner_join' operation in line 64, to make rows unique
            return(res)}
          ) %>% t  #to transpose the matrix
    })
    
 
    
    # assign a letter to each allele in each sample. This will be used later to plot the alleles
    NGSdata.indelranges <- mapply(function(a,b) {
                                    a <- as.data.frame(a)
                                    a$indelID <- letters[1:nrow(a)]
                                    a$sample <- b
                                    rownames(a) <- NULL
                                    return(a)},
                                  NGSdata.indelranges, names(NGSdata.indelranges),
                                 SIMPLIFY = F
                                 )
    
    # indel ranges
    NGSdata.indelranges.ranges <-  do.call(rbind, NGSdata.indelranges) %>%
      mutate(start = Position.IRanges, end = Position.IRanges + indel.IRanges.width)
    
    return(NGSdata.indelranges.ranges)
}

# calculate data for exploratory plots
fun.NGS_exploratory <- function(NGSdata) {
  # NGSdata: data read from the database, produced by fun.NGS_readDB(db)
  
  # This function assesses the samples to decide whether their mutations are homozygous (homo), one allele 'mutant' and the other a different mutant
  #' (muthet), or one allele mutant and the other wt (mutwt). This may need to be adapted for cases with >2 alleles.
      fun.calczygosities <- function(a, offtargets) {
        
        #calculate zygosities
        zygosities <- data.frame(
          'homo' = a$allelefraction >0.85 & nrow(a) == 1,
          'muthet' = a$allelefraction > 0.35 & a$allelefraction <0.85 & nrow(a) == 2,
          'mutwt' = a$allelefraction > 0.35 & a$allelefraction <0.85 & nrow(a) == 1,
          'has.offtargets' = offtargets)
        
        # consolidate zygosities in a single row
        homo <- all(zygosities$homo)
        muthet <- as.logical(sum(zygosities$muthet))   #if it has three rows, two T and one F (e.g. it's a muthet with an off target), it will show as muthet
        mutwt <- as.logical(sum(zygosities$mutwt))
        has.offtargets <- all(zygosities$has.offtargets)
        zygosities <- data.frame(homo, muthet, mutwt, has.offtargets)
        
        # resolve zygosity in a single column. Zygosity will be 'uncertain' if e.g. allelefraction < 0.35
        zygosity.type <- c('homo', 'muthet', 'mutwt')[as.logical(zygosities[1:3])]
        if(length(zygosity.type) == 0) zygosities$zygosity <- 'uncertain' else zygosities$zygosity <- zygosity.type
    
        return(zygosities)
      }
      
  # this function consolidates zygosities in one row per gene    
    fun.calczygosities.checkofftargets <- function(NGSdatalist) {
      ## check for presence of off-targets
      if(length(unique(NGSdatalist$gene)) ==1) { # no off-targets, flag as FALSE    ######## instead of length, use the gene names from target_name. Note that only offtargets that occur at the same time as targets are considered
          zygosities <- fun.calczygosities(NGSdatalist, offtargets = FALSE) %>%
            mutate('gene' = unique(NGSdatalist$gene))
      } else {                         # off-targets
          zygosities <- by(NGSdatalist, INDICES = NGSdatalist$gene, FUN = fun.calczygosities, offtargets = TRUE)
          zygosities <- mapply(function(zygosities.names, zygosities.data) {
                                        zygosities.data$gene <- zygosities.names
                                        return(zygosities.data)},
                               names(zygosities), zygosities, SIMPLIFY = F) %>%
            do.call(rbind,.)
      }
      
      return(zygosities)
    }
    
  # split NGSdata into a list by 'sample' and 'type' (SNV or INDEL)
    NGSdata_list <- split(NGSdata, list(NGSdata$sample, NGSdata$type, NGSdata$guide))
    
  # filter out combinations of factors 'sample' and 'type' that give no rows in the list
    filter_sample.type <- lapply(NGSdata_list, function(a) nrow(a) != 0) %>% do.call(c,.)
    NGSdata_list <- NGSdata_list[filter_sample.type]
  
  # get zygosities
    NGSdata_zygosities <- lapply(NGSdata_list, fun.calczygosities.checkofftargets) %>%
      do.call(rbind,.) %>%
      mutate('ids' = rownames(.))
    
  return(NGSdata_zygosities)
}

# function to calculate percentages of element counts in a variable, used in the plots
fun.percent <- function(a) {
  table(a) %>%
  prop.table(margin = 1) %>%
  as.data.frame() %>%
  mutate(Freq = round(100*Freq, 2))
}

# -----------------------------
# PLOTS

# indel positional ranges
fun.NGS_plotindelranges <- function(NGSdatafull) {

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
  

  # split data into guides (because guides have different locations and they mess with xaxis ranges)
  fun.shortnumber <- function(a) { #this function shortens coordinate numbers into a shortened, plottable version (e.g. 45678987 to 987) 
    a %>%
      as.character() %>%
      substr(nchar(.)-2, nchar(.)) %>%
      as.integer()
  }
  
  NGSdatafull <- subset(NGSdatafull, has.offtargets == F) %>%   #do not plot samples with offtargets
    mutate(short.start = fun.shortnumber(start), short.end = fun.shortnumber(end), fun.shortnumber(start))
  data.perguide <- split(NGSdatafull, NGSdatafull$guide, drop = T)
  
  fun.plot <- function(a) {
    #create plot with no synonymous mutations
    p <-  ggplot() +
      geom_rect(data = subset(a, type == 'INDEL'),
                mapping = aes(xmin = short.start, xmax = short.end, ymin = 0, ymax = allelefraction), fill = 'green4', color = 'grey20') +
      geom_point(data = subset(a, type == 'SNV'),
                mapping = aes(x = short.start , y = allelefraction), color = 'black') +
      geom_vline(xintercept = fun.shortnumber(unique(a$guide_location)), color = 'red', linetype = "dotted", size = 1) +
      xlim((min(a$short.start) - 5), (max(a$short.end) + 5)) +
      xlab('coordinate') +
      theme_bw() +
      theme(strip.text.y = element_text(angle = 0),
            strip.background = element_rect(colour = 'white')) +
      facet_grid(sample ~ guide, scales = 'free_x')
    
    return(p)
  }
  
  return(lapply(data.perguide, fun.plot) %>%
    do.call(grid.arrange,.))
    
}

# Type of mutation (INDEL vs SNV, per guide)
fun.NGS_plotmutations <- function(NGSdatafull) {
  
  mutations <- NGSdatafull %>%
    select(guide, type) %>%
    fun.percent()
  
  plot_ly(mutations, type = 'scatter', mode = 'markers', x = ~type, y = ~Freq, color = ~guide) %>%
    layout(xaxis = list(title = 'Type of mutation'), yaxis = list(title = '%'))
  
}

# Type of variant (frameshift, in-frame deletion, insertion, etc).
fun.NGS_plotvariants <- function(NGSdatafull) {
  
  #shorten names of variants, so they fit better in the x axis
  variants <- NGSdatafull %>%
    select(guide, consequence) %>%
    mutate(consequence = c('frameshift', 'synonymous', 'inframe.del', 'intron.noncoding', 'stop.frameshift')[
    match(NGSdatafull$consequence, c('frameshift', 'synonymous','inframe_deletion', 'intron,non_coding_transcript', 'stop_gained,frameshift'))]) %>%
    fun.percent()
    
  plot_ly(variants, type = 'scatter', mode = 'markers', x = ~consequence, y = ~Freq, color = ~guide) %>%
  layout(xaxis = list(title = 'Type of variant'),
         yaxis = list(title = '%'),
         margin = list(l = 150))
}

# Zygosity (homo, muthet, mutwt). 
  # Note: a muthet has two mutated alleles, each with a different mutation. A mutwt has one mutant allele and one wt allele

fun.NGS_plotzygosity <- function(NGSdatafull) {
  
  zygosity <- NGSdatafull %>%
    select(guide, homo, mutwt, muthet, has.offtargets) %>%
    group_by(guide) %>%
    summarise('total' = n(), 'homo' = sum(homo), 'mutwt' = sum(mutwt), 'muthet' = sum(muthet), 'offtargets' = sum(has.offtargets)) %>%
    mutate_at(vars(homo, muthet, mutwt, offtargets), funs(round(.*100/total, 2))) %>%
    select(-total) %>%
    melt(id.vars = 'guide') %>%
    rename(zygosity = variable, Freq = value)
  
  plot_ly(zygosity, type = 'scatter', mode = 'markers', x = ~zygosity, y = ~Freq, color = ~guide) %>%
    layout(xaxis = list(title = 'Zygosity'), yaxis = list(title = '%'))
}

# Distance of indel to cut site
fun.NGS_plotdistance <- function(NGSdatafull) {
  
  indellength <- NGSdatafull %>%
    select(guide, indellength) %>%
    fun.percent()
  
  fun.plot <- function(a) plot_ly(a, type = 'bar', x = ~indellength, y = ~Freq, color = ~guide) %>%
    layout(xaxis = list(title = 'Indel length (nt)'), yaxis = list(title = '%'))
  
  split(indellength, indellength$guide) %>%
    lapply(fun.plot) %>%
    subplot(shareY = T)
}

# Allele types (e.g. GAA/A)
fun.NGS_plotalleles <- function(NGSdatafull) {
  
  fun.indeltype <- function(a) {
    if(a > 0) b <- 'insertion'
    if(a == 0) b <- 'SNV'
    if(a < 0) b <- 'deletion'
    return(b)
  }
  
  alleles <- NGSdatafull %>%
    select(guide, alleles, indellength, type) %>%
    mutate(indellength = coalesce(indellength, 0L)) %>%  #substitute NA's by 0's
    mutate(indeltype = sapply(indellength, fun.indeltype)) %>%
    mutate(guide = as.factor(guide), alleles = as.factor(alleles))
 
  frequencies <- table(alleles$guide, alleles$alleles) %>%
    prop.table(margin = 1) %>%
    as.data.frame() %>%
    mutate(Freq = round(100*Freq, 2)) %>%
    rename(guide = Var1, alleles = Var2)
  
  alleles <- inner_join(alleles, frequencies, by = c('guide', 'alleles')) %>%
    arrange(indellength) %>%
    mutate(alleles = factor(alleles, levels = unique(alleles)), indeltype = factor(indeltype))
  
  plot_ly(alleles, type = 'scatter', mode = 'markers', x = ~Freq, y = ~alleles, color = ~guide) %>%
    layout(xaxis = list(title = '%'),
           yaxis = list(title = 'Alleles (reference/mutant)'),
           margin = list(l = 250))
}