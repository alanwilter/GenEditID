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
  
  
  
  
  
  
  
  
########Read incuCell data, calculate ratio 800/700 and tabulate:
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

#

# ggplotly(ggplot(subset(incu2, Plate == 'plate5' & Well == 'E3'), aes(x = Elapsed, y = Confluence..., group = Well:Plate:Layout)) +
#   geom_line(aes(colour = Layout)) +
#   theme_bw()
# )

#####PLOTS


##geom_dotplots
        # plot.icw <- ggplot(icw.ratio.mapped, aes(x = Plate, y = ratio800to700)) +
        #   geom_violin() +
        #   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, aes(fill=Layout)) +
        #   theme_bw()
        # 
        # plot.icw2 <- ggplot(icw2, aes(x = Plate, y = ratio800to700)) +
        #   geom_violin() +
        #   geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, aes(fill=Layout)) +
        #   theme_bw()



