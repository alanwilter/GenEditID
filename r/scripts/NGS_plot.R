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

fun.NGS_readDB <- function(db) {

    sql <- 
        'select
             vr.variant_type as Type,
             slc.sequencing_sample_name as Sample,
             slc.sequencing_barcode as Barcode,
             vr.chromosome as Chromosome,
             vr.position as Position,
             vr.allele_fraction as AlleleFraction,
             vr.depth as Depth,
             vr.amplicon as Amplicon,
             vr.indel_length as IndelLength
         from
             sequencing_library_content slc
             inner join variant_result vr on vr.sequencing_library_content_id = slc.id'
    
    dbGetQuery(db$con, sql)
}


#-----------------------------
# CALCULATIONS

# calculate indel positional ranges for fun.NGS_plotindelranges
fun.NGS_indelranges <- function(NGSdata) {
  # NGSdata: data read from the database

    # ---------Ruben update this once the database table is clear  
    # take results from Haplotype.caller only
    #data <-  droplevels(NGSdata[grep('H.', NGSdata$type), ]) # split is done by levels, so make sure that unused ones are dropped
    # ---------
    
    databy <- split(NGSdata, f = NGSdata$sample)
    
    # create indel position ranges (position, indel length), shifting positions for deletions
    # (e.g. an indel in position 8 and width = -2 is converted into indel in postion 6 and width = 2) for IRanges to work
    # Note that in the second line, column 5 is the position and column 9 is the indel length (b[c(5,9)]).
    databy.haploranges <- lapply(databy, function(b) {
      dat <- apply(b[c(5,9)], 1, FUN = function(a) {    #note that `apply` works on matrices, so it must receive just numeric data, otherwise it will convert 'b' into character.
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
    data <- NGSdata[order(NGSdata$sample, NGSdata$position),]
    databy.haploranges.rangesfull <- data.frame(data, databy.haploranges.ranges)
    
    # delete columns 'Sample.1' and 'Original.position'
    databy.haploranges.rangesfull <- databy.haploranges.rangesfull[, -grep('Sample.1|Original.position', colnames(databy.haploranges.rangesfull)),]

    return(databy.haploranges.rangesfull)
}

# calculate data for exploratory plots
fun.NGS_exploratory <- function(db) {
  # NGSdata: data read from the database
  
    sql <- 
        'select
             p.geid as Plate,
             concat(w.row, w.column) as Well,
             g.name as Guide,
             slc.sequencing_sample_name as Sample,
             slc.sequencing_barcode as Barcode,
             vr.consequence as Consequence,
             vr.gene_id as SymbolGene,
             vr."cDNA_effect" as CDNAEffect,
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
             inner join guide g on gwca.guide_id = g.id'
    
    data <- dbGetQuery(db$con, sql)
    
    data <- data[data$allelefraction > 0.15, c('plate', 'well', 'guide', 'alleles', 'indellength',
                                               'allelefraction', 'chromosome', 'position',
                                               'consequence', 'symbolgene')]
    data$content <- paste0('MCF7 clone3 ', data$guide)
    data <- data[order(data$plate, data$well), ]
  
    sql <- 
        'select
             p.geid as Plate,
             concat(w.row, w.column) as Well,
             g.name as Layout,
             slc.sequencing_sample_name as Sample,
             vr.consequence as Variant,
             vr.gene_id as Gene,
             vr.variant_type as Type,
             vr.allele_fraction as AlleleFraction
         from
             well w
             left join sequencing_library_content slc on slc.well_id = w.id
             inner join variant_result vr on vr.sequencing_library_content_id = slc.id
             inner join experiment_layout el on w.experiment_layout_id = el.id
             left join plate p on p.experiment_layout_id = el.id
             inner join well_content wc on w.well_content_id = wc.id
             left join guide_well_content_association gwca on gwca.well_content_id = wc.id
             inner join guide g on gwca.guide_id = g.id'
 
    NGSdata <- dbGetQuery(db$con, sql)
  
    NGSdata <- NGSdata %>%
        mutate('content' = as.factor(paste0('MCF7 clone3 ', layout)),
               'plate' = as.factor(plate)) %>%
        mutate('plate' = gsub('plate', '', plate)) %>%
        droplevels %>%
        arrange(desc(plate, well))
  
  NGSdata.cells <- NGSdata[grepl('-C', NGSdata$sample),] # extracted cells only
  NGSdata.gDNA <- NGSdata[grepl('-G', NGSdata$sample),]  # extracted gDNA only
  
  fun.by <- function(x) {NGS <- by(x[,], INDICES = list(x$plate, x$well, x$type), FUN = as.data.frame)
  NGS <- NGS[!do.call(c, lapply(NGS, is.null))]
  NGS <- lapply(NGS, function(a) data.frame(
    'Plate' = a$plate,
    'Well' = gsub('P.', '', a$well),
    'Content' = paste0('MCF7 clone3 ', a$layout),
    'Type' = a$type,  #SNV or INDEL
    'NROW' = nrow(a), 
    'homo' = a$allelefraction >0.85,
    'muthet' = a$allelefraction > 0.35 & a$allelefraction <0.85 & nrow(a) > 1,
    'mutwt' = a$allelefraction > 0.35 & a$allelefraction <0.85 & nrow(a) == 1,
    'has.offtargets' = length(unique(a$gene)) > 1)
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
  
  fun.by2 <- function(x) {NGS <- by(x[,], INDICES = list(x$plate, x$well), FUN = as.data.frame)
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
  
  
g.plot <- ggplot(subset(a, indel.IRanges.width > 0), aes(xmin = start, xmax = end, ymin = 0, ymax = allelefraction, color = indelID)) +
  geom_point(data = subset(a, indel.IRanges.width == 0 & chromosome == 'chr10'), aes(x = position, y = allelefraction), color = 'black') +
  geom_rect(mapping = aes(fill = type)) +
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
