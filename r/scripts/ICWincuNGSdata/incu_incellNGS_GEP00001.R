library('reshape2')
library('plotly')
library('grofit')
library('pheatmap')

###Functions
func.reorder <- function(a) do.call(rbind, as.list(a))
func.extractplate <- function(a) {
  plate <- sapply(a, FUN = function(b) strsplit(b, '_|\\.'))
  plate <- lapply(plate, function(b) tolower(b[grep('plate', b, ignore.case = T)]))
  plate <- do.call(c, plate)
}

###Values
wells <- func.reorder(outer(LETTERS[1:8], 1:12, paste0))

##########BODY##########
###Read plates layout
plate.layout <- as.matrix(read.delim('PlatesLayout_270117.csv', header = F))
pln <- tolower(unique(as.character(plate.layout[,15]))) #pln stands for 'plate layout number'
pln <- pln[!pln == '']
seq.pln <- 0:(length(pln)-1)
plate.layout <- lapply(seq.pln, function(a) plate.layout[c((2 + (11*a)) : (9 + (11*a))), 2:13])
plate.layout <- do.call(rbind,
                        mapply(function(a, b) data.frame('Plate' = a, 'Well' = wells, 'Layout' = func.reorder(b)),
                               pln, plate.layout, SIMPLIFY = F)
)

###Read incellWestern data, calculate ratio 800/700 and tabulate:
#Get file names
files.icw <- list.files(pattern='ICW', ignore.case=T)
#Get plate numbers from file names (add warning in shiny if no 'plate' is found on the file name)
files.icw.plate <- func.extractplate(files.icw)
#Get data and process it
icw700 <- lapply(files.icw, function(a) as.matrix(read.delim(a, skip=24, nrows=8, header=F)))
icw800 <- lapply(files.icw, function(a) as.matrix(read.delim(a, skip=34, nrows=8, header=F)))
icw.ratio <- mapply(function(a,b) a/b, icw800, icw700, SIMPLIFY = F) #this makes a column matrix, ordered as A1, B1, C1...A2, B2, C2...

icw.ratio <- do.call(rbind,
                     mapply(function(a1, a2, a3, a4, a5) data.frame('PlateID.ICW' = a1, 
                                                            'Plate' = a2,
                                                            'Well' = wells,
                                                            'ratio800to700' =func.reorder(a3), 
                                                            'read700nm' = func.reorder(a4),
                                                            'read800nm' = func.reorder(a5)),
                          files.icw, files.icw.plate, icw.ratio, icw700, icw800, SIMPLIFY = F)
             )
icw.ratio.mapped <- merge(plate.layout, icw.ratio, by = c('Plate', 'Well'), all = T, sort = F)

#Which ratio values are just background?
 #create a variable that shows empty vs filled wells
icw2 <- cbind(icw.ratio.mapped, 'mp' = icw.ratio.mapped$Layout == '-')

 #plots before cleaning empty wells' background
  ggplotly(ggplot(icw2, aes(x = Plate, y = ratio800to700)) +
             geom_violin() +
             geom_point(aes(group=Well, colour=mp)) +
             theme_bw()
  )
  
  ggplotly(ggplot(icw2, aes(x = Plate, y = read800nm)) +
             geom_violin() +
             geom_point(aes(group=Well, colour=mp)) +
             theme_bw()
  )
  
  gp <- ggplotly(ggplot(icw2, aes(x = Plate, y = read700nm)) +
             geom_violin() +
             geom_point(aes(group=Well, colour=mp)) +
             theme_bw()
  )
  
 #get rid of all empty wells that have read700nm values that are not outliers (i.e. to get rid of wells with cells that failed to grow)
  icw2.empty <- icw2[icw2$mp == T, 'read700nm']
  icw2.empty2 <- mean(icw2.empty) + 2*sd(icw2.empty)
  icw2.empty2 <- icw2.empty[icw2.empty > icw2.empty2]
  icw2 <- icw2[! (icw2$mp == T & icw2$read700nm < min(icw2.empty2)),]

  # extract values from MCF7_ctl to see if there is a positional effect
  icw2.MCF7 <- subset(icw2, Plate == 'plate6')
  # icw2.MCF7Layout <- matrix(icw2.MCF7[, 3], nrow = 8, ncol = 12)
  # icw2.MCF7values <- matrix(icw2.MCF7[, 5], nrow = 8, ncol = 12)
  icw2.MCF7$WellA <- substr(icw2.MCF7$Well, 0, 1)
  icw2.MCF7$Welln <- substr(icw2.MCF7$Well, 2, 3)
  
  ###################################I need to do this with the full plate before trimming
  ggplot(icw2.MCF7, aes(x = as.numeric(Welln), y = WellA)) + 
    geom_tile(aes(fill = read700nm), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(label=Layout), size = 3) +
    scale_x_continuous('Well', breaks = 1:12) +
    scale_y_discrete(limits = rev(levels(as.factor(icw2.MCF7$WellA))))
  
  
#Plots with background-corrected values
  ggplotly(ggplot(icw2, aes(x = Plate, y = read700nm)) +
             geom_violin() +
             geom_point(aes(group=Well, colour=Layout)) +
             #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, aes(fill=Layout)) +
             theme_bw()
  )
  
  ggplotly(ggplot(icw2, aes(x = Plate, y = read800nm)) +
             geom_violin() +
             geom_point(aes(group=Well, colour=Layout)) +
             # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, aes(fill=Layout)) +
             theme_bw()
  )
  
  ggplotly(ggplot(icw2, aes(x = Plate, y = ratio800to700)) +
             geom_violin() +
             geom_point(aes(group=Well, colour=Layout)) +
             # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, aes(fill=Layout)) +
             theme_bw()
  )

  ggplotly(ggplot(icw2, aes(x = Layout, y = ratio800to700)) +
             geom_violin() +
             geom_point(aes(colour=Plate)) +
             # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, aes(fill=Layout)) +
             theme_bw()
  )
  
#subset(icw2, icw2$Plate == 'plate1' & icw2$Well %in% c('D10', 'F4', 'E6','C10'))
#subset(icw2, icw2$Plate == 'plate1' & icw2$Well %in% c('C10', 'B3', 'B10', 'F10'))
#subset(icw2, icw2$Plate == 'plate4' & icw2$Well %in% c('E5'))
  #Show potential Knock-downs and KO's (i.e. cells that have a ratio800/700 smaller than the MCF7 control)
   #get rid of potential outliers in the MCF7 control, ie ratio values < mean + 2*sd, and set the threshold of homozygous
   #'expression at the minimum of the MCF7 control values that are left.
  
  
  
  
  
  
  
  
########Read incucyte data, calculate ratio 800/700 and tabulate:
#Get file names
files.incu <- list.files(pattern='luc-straw', ignore.case=T)
#Get plate numbers from file names (add warning in shiny if no 'plate' is found on the file name)
files.incu.plate <- func.extractplate(files.incu)
incu <- lapply(files.incu, read.delim, skip = 7)
incu <- mapply(function(a1, a2, a3) data.frame('PlateID.IC' = a1, 'Plate' = a2,
                                         melt(a3, id.vars=1:2, measure.vars = 3:ncol(a3), variable.name = 'Well', value.name = 'Confluence(%)')),
               files.incu, files.incu.plate, incu, SIMPLIFY = F)

incu <- do.call(rbind, incu)
#map Layout
incu <- merge(plate.layout, incu, by = c('Plate', 'Well'), all = T, sort = F)

# Trim wells with growth increase (time 93 - time 0) < 10%
incu.growthincrease <- incu[incu$Elapsed == 93, 'Confluence...'] - incu[incu$Elapsed == 0, 'Confluence...']
incu.growthincrease <- cbind(incu[incu$Elapsed == 0, c('Plate', 'Well')],
                             'GrowthIncrease' = incu.growthincrease, 'GrowthincreaseMoreThan10' = incu.growthincrease > 12)
      #it could be worth it to make the growthincrease threshold manually selectable (so a user can say 'I prefer faster growers, etc')

incu.growthincrease <- incu.growthincrease[incu.growthincrease$GrowthincreaseMoreThan10 == T,]

incu.index <- paste0(incu$Plate, incu$Well) %in% paste0(incu.growthincrease$Plate, incu.growthincrease$Well)
incu2 <- incu[incu.index,]




#Convert incu2 format to matrix to use with the grofit package
incu2.data <- dcast(incu2, Plate + Well + Layout ~ Elapsed, value.var = 'Confluence...') #NOTE: the way grofit() is written, the .data must have 3 more 
incu2.time <- matrix(data = incu2$Elapsed, ncol = length(unique(incu2$Elapsed)), nrow = nrow(incu2.data), byrow = T)
incu.grofit <- grofit(incu2.time, incu2.data, control = grofit.control(interactive=F))
incu.grofitsum <- summary(incu.grofit$gcFit)
incu.grofitsum <- incu.grofitsum[,c(1:9, 13)]
colnames(incu.grofitsum)[c(1:3, 9, 10)] <- c('Plate', 'Well', 'Layout', 'mu', 'stdmu') #mu is the maximum slope, and stdmu its std deviation

#plot clone growth rates
ggplotly(ggplot(incu.grofitsum, aes(x = Layout, y=mu, group = Layout:Plate:Well)) +
  geom_point(position = position_dodge(0.2), aes(colour = Layout)) +
  geom_errorbar(aes(ymin = mu - stdmu, ymax= mu + stdmu, colour = Layout), position = position_dodge(0.2), size = 0.4) +
  theme_bw() +
  ggtitle('Clone growth rate') +
  xlab('Clone') +
  ylab('maximum growth slope') +
  theme(
    #axis.title.x = element_text(face="bold", colour="#990000", size=20),
    axis.text.x  = element_text(angle=15, vjust=0.5, size=10))
)

subset(incu.grofitsum, Well == 'F10' & Plate == 'plate5')

#plot clone growth
ggplotly(ggplot(incu2, aes(x = Elapsed, y = Confluence..., group = Plate:Well:Layout)) +
  geom_line(aes(colour = Layout)) +
  theme_bw()
)

ggplotly(ggplot(subset(incu2, Plate == 'plate3'& Well == 'E4'), aes(x = Elapsed, y = Confluence..., group = Plate:Well:Layout)) +
           geom_line(aes(colour = Layout)) +
           theme_bw()
)
#


#merge ICW and incu data in a single table for plotting:
ICWincu <- merge(incu.grofitsum, icw2, by = 1:3)
#this is a list with the KO clones selected by Rasmus
KOlist <- data.frame('Plate' = paste0('plate', c(1,1,1,2,2,3,3,3,4,4)),
                     'Well' = c('C5', 'C9', 'E2', 'D9', 'F8', 'B9', 'C3', 'D9', 'B6', 'B5'))

ICWincu$mer <- paste0(ICWincu$Plate, ICWincu$Well)
KOICWincu <- subset(ICWincu, mer %in% paste0(KOlist[,1], KOlist[,2]))

subset(ICWincu, Plate %in% KOlist[,1] & Well %in% KOlist[,2])

ggplotly(ggplot(ICWincu, aes(x = ratio800to700, y=mu, group = Layout:Plate:Well)) +
           geom_point(aes(colour = Layout)) +
           #geom_errorbar(aes(ymin = mu - stdmu, ymax= mu + stdmu), position = position_dodge(0.2), size = 0.4) +
            #NOTE!! position_dodge messes up with the hover events and gives wrong locations, don't use!
           geom_point(data = KOICWincu, colour= 'black') +
           theme_bw() +
           ggtitle('Clone growth rate') +
           xlab('protein abundance') +
           ylab('maximum growth slope') +
           theme(
             #axis.title.x = element_text(face="bold", colour="#990000", size=20),
             axis.text.x  = element_text(angle=15, vjust=0.5, size=10))
)



###Add in NGS data
    snv <- read.delim("SLX-13775_variantsSNVs.csv")
    snv$Type <- 'SNV'
    gDNAsnv <- sort(unique(as.character(snv[grep('*-G', snv[,1]), 1])))
    
    indel <- read.delim("SLX-13775_variantsINDELS.csv")
    indel$Type <- 'INDEL'
    gDNAindel <- sort(unique(as.character(indel[grep('*-G', indel[,1]),1])))
    
    gDNAcontrols <- as.character(read.delim("gDNAlines.csv")[,1])
    
    #A few checks to see if all gDNA controls are represented in the variant calling data
    data <- rbind(snv, indel)
    setdiff(gDNAsnv, gDNAindel)
    setdiff(union(gDNAsnv, gDNAindel), gDNAcontrols)
    un <-union(gDNAsnv, gDNAindel)
    setdiff(un, gDNAcontrols)
    
    
    #Add guide notations to 'data'
    data$guide <- substr(as.character(data[,1]), 4, 5)
    data$guide <- c('STAT3.3','STAT3.4','STAT3.1','STAT3.2')[match(data$guide,c('P1', 'P2', 'P3', 'P4'))] 
    data$Plate <- substr(data$Sample, 4,5)
    data$Plate <- gsub('P', 'plate', data$Plate)
    data$Plate <- as.factor(data$Plate)
    data$Well <- substr(data$Sample, 6,7)
    data$Well <- as.factor(data$Well)
    dataNGS <- data[,c(36, 37, 35, 1, 3, 4, 34, 15)]
    colnames(dataNGS) <- c('Plate', 'Well', 'Layout', 'Sample', 'Variant', 'Gene', 'Type', 'AlleleFraction')
    dataNGS$Type <- as.factor(dataNGS$Type)
    
    #Filter out rows with a very low AlleleFraction (i.e. false positive SNVs and indels due to PCR or sequencing errors,
    #' because normally we expect an allele fraction of 1/2 for heterozygotes and 1 for homozygotes in biallelic cases.
    #' Do not be too stringent here, because e.g. in triallelic cases we would expect a fraction of 1/3 or above for ranges
    #' of heterozygotes to homozygotes)
    dataNGS <- dataNGS[dataNGS$AlleleFraction > 0.15,]
    
#Merge ICW, IncuCyte and NGS data
  ICWincu <- ICWincu[,c(1:3, 9, 10, 12)]
  ICWincu <- ICWincu[grepl('Clone|STAT',ICWincu$Layout),]
  ICWincuNGS <- merge(ICWincu, dataNGS, by = c('Plate', 'Well', 'Layout'), all.x = T)
  ICWincuNGS$plusSD <- ICWincuNGS$mu + ICWincuNGS$stdmu
  ICWincuNGS$minusSD <- ICWincuNGS$mu - ICWincuNGS$stdmu

  
# It would be nice to have this plot linked to the NGS table, so when hovering on a point or clicking on it the table comes up
#' and we can see the more detailed data for the clone (e.g. frameshifts, allele fractions, etc.).
#' Coding for allele fractions in a visual way would be nice too (e.g. AlleleFraction around 0.5, heterozygote; fraction around 1, homozygote);
#' but we ran into the problem of multiallelic cases, where these automatic calls would not work.
#' I have exported this plot as ICWincuNGS.html.
#' Note that the ICW and Incu data represents more clones than the NGS data, so in a more refined plot below I show the ICW and incu data
#' for only the clones screened by NGS.
ggplotly(ggplot(ICWincuNGS, aes(x = ratio800to700, y=mu, group = Plate:Well:Layout)) +
           geom_point(aes(shape = Layout), size = 1) +       #each guide gets a shape. Probably it is too much cluttering. For an exploratory plot
                                                             #' this is fine
           geom_errorbar(aes(ymin = minusSD, ymax= plusSD), size = 0.4) +
           geom_point(data = KOICWincu, colour= 'green') +   #library sample names ending in -G, i.e. clones gDNA was extracted from.
           geom_point(data = ICWincuNGS[!is.na(ICWincuNGS$Type),], aes(color = Type)) + # indels and snvs. Note overlaps in clones, 
                                                              #'so we need to check IGV data to see which ones are real (probably the SNVs
                                                              #'are included in the indels)
           theme_bw() +
           ggtitle('Clone growth rate') +
           xlab('protein abundance') +
           ylab('maximum growth slope') + #calculated with the grofit package. The vertical bars denote range of confidence for the slope.
                                          #'hopefully the curve smoothing for growth data will get better values.
           theme(
             axis.text.x  = element_text(angle=15, vjust=0.5, size=10))
)

#Plot for only NGS screened clones
#' It shows that the efficiency of editing is very high! Most of the indels generate low-protein mutants

NGSscreenedclones <- read.delim('NGSscreenedClones.csv', stringsAsFactors = F, header = F) #These are the clones that went for NGS
NGSscreenedclones <- melt(NGSscreenedclones, id.vars = NULL)[,2]
NGSscreenedclones <- cbind(substr(NGSscreenedclones, 1, 6), substr(NGSscreenedclones, 7, 8))
ICWincuNGSonly <- merge(NGSscreenedclones, ICWincuNGS, by = 1:2, all.y = F, sort = F)
colnames(ICWincuNGSonly)[1:2] <- c('Plate', 'Well')
ICWincuNGSonly <- rbind(ICWincuNGSonly, subset(ICWincuNGS, Layout == 'Clone3')) #ad Clone3 (the WT control) for protein and growth reference

ICWincuNGSonly$AlleleFractionType <- 'Het'
ICWincuNGSonly[!is.na(ICWincuNGSonly$Gene) & ICWincuNGSonly$AlleleFraction > 0.8, 'AlleleFractionType'] <- 'Homo'
ICWincuNGSonly[is.na(ICWincuNGSonly$Gene), 'AlleleFractionType'] <- 'Control'
ICWincuNGSonlyGenes <- ICWincuNGSonly[!is.na(ICWincuNGSonly$Gene),]
             
#I still need to get colour mapping right, it's not easy to interpret unless you keep on clicking on the legend.
ggplotly(ggplot(ICWincuNGSonly, aes(x = ratio800to700, y=mu, group = Plate:Well:Layout)) +
           geom_point(size = 1) +
           geom_errorbar(aes(ymin = minusSD, ymax= plusSD), size = 0.4) +
           geom_point(data = KOICWincu, colour= 'green') +   #library sample names ending in -G, i.e. clones gDNA was extracted from.
           geom_point(data = ICWincuNGSonlyGenes, aes(colour = Type)) + # indels and snvs. Note overlaps in clones, 
           geom_point(data = ICWincuNGSonlyGenes, aes(colour = AlleleFractionType)) + 
               #'so we need to check IGV data to see which ones are real (probably the SNVs
           #'are included in the indels)
           theme_bw() +
           ggtitle('Clone growth rate') +
           xlab('protein abundance') +
           ylab('maximum growth slope') + #calculated with the grofit package. The vertical bars denote range of confidence for the slope.
           #'hopefully the curve smoothing for growth data will get better values.
           theme(
             axis.text.x  = element_text(angle=15, vjust=0.5, size=10))
)
