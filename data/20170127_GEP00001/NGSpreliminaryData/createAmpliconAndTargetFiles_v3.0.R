# This script is modified to fit to new format 
# data 05-04-2017

library( dplyr)



makeTargetsUniq <- function( tab )
{
  tab <- mutate(tab, IDs = paste( chromosome, start, end, sep=':'))
  uniqueTab <- tab[!duplicated(tab$IDs),]
  
  for (EachDupId in unique(tab$IDs[duplicated(tab$IDs)]) ){
    #cat( EachDupId, "\n")
    newName <- paste(as.vector(filter(tab, IDs == EachDupId)$ampliconName), collapse = ',' )
    uniqueTab[uniqueTab$IDs == EachDupId, 'ampliconName'] <- newName
  }
  
  uniqueTab <- mutate(uniqueTab, chromosome= paste( 'chr', chromosome, sep=''))
  uniqueTab <- arrange(uniqueTab, chromosome, start)
  uniqueTab <- select(uniqueTab, -IDs)
  return(uniqueTab)
  
}




tab <- read.csv( file='20170118_GEP00001_format_05_04_17.csv', stringsAsFactors = F)
tab <- select(tab, -description)

strand <- '-'
chromosome <- c()
ampliconStart <- c()
ampliconEnd <- c()
targetStart <- c()
targetEnd <- c()
ampliconName <- c()


for( amplicon in unique(tab$amplicon))
{
  cat( amplicon, "\n")
  
  # get amplicon start
  aStart <- filter(tab, amplicon_name == amplicon, primer_strand == 'forward')$primer_start
  ampliconStart <- c( ampliconStart, aStart)
  
  # get amplicon end
  aEnd <- filter(tab, amplicon_name == amplicon, primer_strand == 'reverse')$primer_end
  ampliconEnd <- c( ampliconEnd, aEnd)
  
  # get target start
  tStart <- filter(tab, amplicon_name == amplicon, primer_strand == 'forward')$primer_end + 1
  targetStart <- c(targetStart, tStart)
  # get target end
  tEnd <- filter(tab, amplicon_name == amplicon, primer_strand == 'reverse')$primer_start - 1
  targetEnd <- c(targetEnd, tEnd)
  
  # get chromosome info
  chr <- filter(tab, amplicon_name == amplicon, primer_strand == 'forward')$chromosome
  chromosome <- c(chromosome, chr)
  
  # get amplicon name
  aName <- filter(tab, amplicon_name == amplicon, primer_strand == 'forward')$amplicon_name
  ampliconName <- c(ampliconName, aName)
  
}

ampliconTab <- data.frame( chromosome=chromosome, start=ampliconStart, end=ampliconEnd, 
                           strand=rep(strand, length(chromosome)), ampliconName=ampliconName, stringsAsFactors = F)
ampliconTabUniq <- makeTargetsUniq(ampliconTab)


write.table(ampliconTabUniq, file = 'amplicons_v1.txt', col.names = F, row.names = F, quote = F, sep='\t')

ampliconTabUniq <- mutate(ampliconTabUniq, start=start-1)

ampliconTabUniq <- select(ampliconTabUniq, -strand)

write.table( ampliconTabUniq, file = 'amplicons_v1.bed', col.names = F, row.names = F, quote = F, sep='\t')



# target files
targetTab <- data.frame( chromosome=chromosome, start=targetStart, end=targetEnd, 
                         strand=rep(strand, length(chromosome)), ampliconName=ampliconName, stringsAsFactors = F)
targetTabUniq <- makeTargetsUniq(targetTab)

write.table(targetTabUniq, file = 'target_v1.txt', col.names = F, row.names = F, quote = F, sep='\t')

targetTabUniq <- mutate(targetTabUniq, start = start -1)

targetTabUniq <- select(targetTabUniq, -strand)
write.table(targetTabUniq, file = 'target_v1.bed', col.names = F, row.names = F, quote = F, sep='\t')

