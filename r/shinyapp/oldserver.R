
#########################################################

###Generic functions
func.reorder <- function(a) do.call(rbind, as.list(a))
func.extractplate <- function(a) {
  plate <- sapply(a, FUN = function(b) strsplit(b, '_|\\.'))
  plate <- lapply(plate, function(b) tolower(b[grep('plate', b, ignore.case = T)]))
  plate <- do.call(c, plate)
}

###Values
wells <- func.reorder(outer(LETTERS[1:8], 1:12, paste0))

###Read plates layout
func.readlayout <- function(a) {
  plate.layout <- as.matrix(read.delim(a$datapath, header = F))
  pln <- tolower(unique(as.character(plate.layout[,15]))) #pln stands for 'plate layout number'
  pln <- pln[!pln == '']
  seq.pln <- 0:(length(pln)-1)
  plate.layout <- lapply(seq.pln, function(a) plate.layout[c((2 + (11*a)) : (9 + (11*a))), 2:13])
  plate.layout <- do.call(rbind,
                          mapply(function(a, b) data.frame('Plate' = a, 'Well' = wells, 'Layout' = func.reorder(b)),
                                 pln, plate.layout, SIMPLIFY = F)
  )
  return(plate.layout)
}

###Read ICW data
func.readICW <- function(a, playout) {
  #a is inFile
  #playout is plate.layout
  #Get file names
  files.icw <- list(a$name)
  #Get plate numbers from file names (add warning in shiny if no 'plate' is found on the file name)
  files.icw.plate <- func.extractplate(files.icw)
  #Get data and process it
  icw700 <- lapply(a$datapath, function(a) as.matrix(read.delim(a, skip=24, nrows=8, header=F)))
  icw800 <- lapply(a$datapath, function(a) as.matrix(read.delim(a, skip=34, nrows=8, header=F)))
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
  icw.ratio.mapped <- merge(playout, icw.ratio, by = c('Plate', 'Well'), all = T, sort = F)
  
  #Which ratio values are just background?
  #create a variable that shows empty vs filled wells
  icw2 <- cbind(icw.ratio.mapped, 'mp' = icw.ratio.mapped$Layout == '-')
  
  #get rid of all empty wells that have read700nm values that are not outliers (i.e. to get rid of wells with cells that failed to grow)
  icw2.empty <- icw2[icw2$mp == T, 'read700nm']
  icw2.empty2 <- mean(icw2.empty) + 2*sd(icw2.empty)
  icw2.empty2 <- icw2.empty[icw2.empty > icw2.empty2]
  icw2 <- icw2[! (icw2$mp == T & icw2$read700nm < min(icw2.empty2)),]
  return(icw2)
}

###Read incuCell data, calculate ratio 800/700 and tabulate:
func.readIncu <- function(a, playout) {
  #a is inFile
  #playout is plate.layout
  
  #Get file names
  files.incu <- a$name
  #Get plate numbers from file names (add warning in shiny if no 'plate' is found on the file name)
  files.incu.plate <- func.extractplate(files.incu)
  incu <- lapply(a$datapath, read.delim, skip = 7)
  inc<<- incu
  incu <- mapply(function(a1, a2, a3) data.frame('PlateID.IC' = a1, 'Plate' = a2,
                                                 melt(a3, id.vars=1:2, measure.vars = 3:ncol(a3), variable.name = 'Well', value.name = 'Confluence(%)')),
                 files.incu, files.incu.plate, incu, SIMPLIFY = F)
  
  incu <- do.call(rbind, incu)
  #map Layout
  incu <- merge(playout, incu, by = c('Plate', 'Well'), all = T, sort = F)
  
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
  inc <<- list(incu2, incu.grofitsum)
  return(list(incu2, incu.grofitsum))
}



func.plotICW <- function(a) { 
  plot1 <- ggplot(a, aes(x = Layout, y = ratio800to700)) +
    geom_violin() +
    geom_point(aes(colour=Layout)) +
    # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, aes(fill=Layout)) +
    theme_bw() + 
    ggtitle('InCellWestern') +
    xlab('Cell line') +
    ylab('Relative protein abundance (ratio 800 to 700 nm)') +
    theme(
      #axis.title.x = element_text(face="bold", colour="#990000", size=20),
      axis.text.x  = element_text(angle=15, vjust=0.5, size=10))
  
  return(plot1)
}

func.plotIncuRate <- function(a) { 
  plot1 <- ggplot(a, aes(x = Layout, y=mu, group = Layout:Plate:Well)) +
    geom_point(position = position_dodge(0.2), aes(colour = Layout)) +
    geom_errorbar(aes(ymin = mu - stdmu, ymax= mu + stdmu, colour = Layout), position = position_dodge(0.2), size = 0.4) +
    theme_bw() +
    ggtitle('Clone growth rate') +
    xlab('Clone') +
    ylab('maximum growth slope') +
    theme(
      #axis.title.x = element_text(face="bold", colour="#990000", size=20),
      axis.text.x  = element_text(angle=15, vjust=0.5, size=10))
  return(plot1)
}

func.plotIncuCurve <- function(a) {
  ggplot(a, aes(x = Elapsed, y = Confluence..., group = Layout:Plate:Well)) +
    geom_line(aes(colour = Layout)) +
    theme_bw() +
    ggtitle('Clone growth curves') +
    xlab('Clone') +
    ylab('Confluence (%)') +
    theme(
      #axis.title.x = element_text(face="bold", colour="#990000", size=20),
      axis.text.x  = element_text(angle=15, vjust=0.5, size=10))
}


shinyServer(function(input, output, session){
  
  #Read files, all ICW and IncuCyte together   
  allfiles <- reactive({
    validate(need(input$fileinput != '', 'Enter one or several csv files'))
    inFile <- input$fileinput
    if (is.null(inFile)) return(NULL)
    #csvfile <- read.csv(inFile$datapath, stringsAsFactors = F, header = F)
    csv <<- inFile
    files.layout <- inFile[grep('layout', inFile$name, ignore.case = T), ]
    validate(need(nrow(files.layout) == 1, 'You need to supply a plate layout tab-delimited file'))
    plate.layout <- func.readlayout(files.layout)
    
    files.ICW <- inFile[grep('ICW', inFile$name, ignore.case = T), ]
    #validate(need(nrow(files.layout) != 0, 'You need to supply a plate layout tab-delimited file'))
    plate.ICW <- func.readICW(files.ICW, plate.layout)
    icw <<- plate.ICW
    
    files.incu <- inFile[grep('MCF7', inFile$name, ignore.case = T), ]
    #validate(need(nrow(files.incu) != 0, 'You need to supply a plate layout tab-delimited file')) 
    plate.incu <- func.readIncu(files.incu, plate.layout)
    inc <<- plate.incu
    
    return(list(plate.ICW, plate.incu))
    
  })
  
  #ICW plot
  ICWplotly <- reactive({
    allfiles()
    plotICW <- allfiles()[[1]]
    plot1 <- func.plotICW(plotICW)
    return(plot1)
  })
  
  # pplot <- reactive({
  #   ppl <- ICWplotly() + plot_ly(source = 'source')
  # })
  
  IncuRateplotly <- reactive({
    allfiles()
    plot1 <- func.plotIncuRate(allfiles()[[2]][[2]])
    return(plot1)
  })
  
  IncuCurveplotly <- reactive({
    allfiles()
    plot1 <- func.plotIncuCurve(allfiles()[[2]][[1]])
  })
  
  #output$tablecsv <- renderTable(allfiles())
  output$ICWplot <- renderPlotly(ICWplotly())
  output$Incucyteplot <- renderPlotly(IncuRateplotly())
  output$IncucyteCurveplot <- renderPlotly(IncuCurveplotly())
  # eventdata <<- isolate(event_data("plotly_hover", source = "source"))
  
  
}
)