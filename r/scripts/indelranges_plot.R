# DESCRIPTION:
# R script to manipulate NGS indel data and create a tile plot with indel positions and allele frequencies per well.
# Currently it does not read from the database. Also we need to modify the hard-added column 'Type', so it is managed from the database only.
# Things that should come or go from/to the database are marked as %DB

library(ggbio)
library(IRanges)
library(dplyr)
library(GenomicRanges)

# Hard-added a column called 'Type' with the variant caller used (Haplotype.caller or Vardict)
data <- read.delim('SLX-13774variants.csv') # I have not put this file on github, but Variant files are in data/GEP00002/ %DB
data <- data[,c(1:4, 9, 10, 16,17, 20, 31)]
colnames(data) <- c('Type', 'Sample', 'Barcode', 'Variant', 'Chromosome', 'Position', 'Allele.fraction', 'Depth', 'Amplicon', 'indel.length')

# take results from Haplotype.caller only
data <-  droplevels(data[grep('H.', data$Type), ]) # split is done by levels, so make sure that unused ones are dropped
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

      #To convert to IRanges object
        # databy.haploranges.ranges <- IRanges(start = databy.haploranges.ranges[,1], width = databy.haploranges.ranges[,2])
        # 
        #   
        # databy.granges <- GRanges(seqnames = Rle('chr10', length(databy.haploranges.ranges)),
        #         ranges = databy.haploranges.ranges,
        #         strand = Rle('+', length(databy.haploranges.ranges)),
        #         sample = Rle(names(databy.haploranges.nalleles), databy.haploranges.nalleles)
        #         )
        # databy.granges <- as.data.frame(databy.granges)
        # databy.granges$end <- databy.granges$end + 1 #(since this is not a GRanges object anymore, I have to adjust widths, so SNVs show as start=1, end=1, and indel as start=1, end=2)

# merge position ranges with original indel data
 # note that sample names are duplicated (because there are several indels per sample), 
 # so `merge` or `inner_join` will fail and generate incorrect (bigger) dataframes.
 # Instead of that, merge by ordering and binding the data.
databy.haploranges.ranges <- databy.haploranges.ranges[order(databy.haploranges.ranges$Sample, databy.haploranges.ranges$Original.position), ]
data <- data[order(data$Sample, data$Position),]
databy.haploranges.rangesfull <- data.frame(data, databy.haploranges.ranges)

# delete columns 'Sample.1' and 'Original.position'
databy.haploranges.rangesfull <- databy.haploranges.rangesfull[, -grep('Sample.1|Original.position', colnames(databy.haploranges.rangesfull)),]

g.plot <- ggplot(subset(databy.haploranges.rangesfull, indel.IRanges.width > 0), aes(xmin = start, xmax = end, ymin = 0, ymax = Allele.fraction, color = indelID)) +
  geom_point(data = subset(databy.haploranges.rangesfull, indel.IRanges.width == 0 & Chromosome == 'chr10'), aes(x = Position, y = Allele.fraction), color = 'black') +
  geom_rect(mapping = aes(fill = Variant)) +
  facet_wrap(~Sample) +
  geom_vline(xintercept=89653802, color = 'steelblue', linetype = "dotted") +
  theme_classic()


# ggplotly(g.plot)

write.table(databy.haploranges.rangesfull, 'PTENdata.txt', sep = '\t', row.names = F) #%DB This should go to the database instead
