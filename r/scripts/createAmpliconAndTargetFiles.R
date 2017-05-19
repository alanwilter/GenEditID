library( 'dplyr')
library( 'RPostgreSQL')

# connect to DB
db <- src_postgres(user="gene", password="gene", host="bioinf-ge001.cri.camres.org", 
                   port=5432, dbname="geneediting" )


primer_amplicon_association <- tbl( db, 'primer_amplicon_association') %>% as.data.frame()

primer <- tbl( db, 'primer') %>% 
  select(id, strand, start, end) %>% 
  rename(primer_id=id,  primer_strand=strand, primer_start=start, primer_end=end) %>% 
  as.data.frame()

amplicon <- tbl(db, 'amplicon') %>% 
  select( id,chromosome, start, end) %>%
  rename( amplicon_id =id, amplicon_start=start, amplicon_end=end) %>% 
  as.data.frame()


target <- tbl(db, 'target') %>% as.data.frame()


  #               <----Amplicon----------------->
  #               ===============================
  #  Far Primer   >>>>>>>          <<<<<<< Rev primer
  #                      <-target->


# get strand information
# Where dose amplicon seq pipeline usese strand information?
# I am not clear, needs to discuss with Matt
# Currently all the coordinates are forward starand based

#strand <- ifelse( target$strand == 'reverse', '-', '+')
strand <- '+'  
genome_id <- target$genome_id


# amplicon coordinates
#  GRC reference genomes chromosome names don't have 'chr'
amplicon <- mutate(amplicon, chromosome = paste( 'chr', chromosome, sep='')) %>% 
  mutate( strand=strand, amplicon_name=paste(genome_id, chromosome, amplicon_start, sep='_' )) %>% 
  arrange(chromosome, amplicon_start) 

write.table(x=select(amplicon, -amplicon_id), file='amplicons_v1.0.txt',
            row.names = FALSE, col.names = FALSE, quote=FALSE, se='\t'
            )

# Target coordinates

# join primer_amplicon_association + amplicon
paa_pl_a <- inner_join(primer_amplicon_association, amplicon, by='amplicon_id')

# paa_pl_a + primer
paa_pl_a_pl_p <- inner_join( primer, paa_pl_a, by='primer_id')



# get target coordinates
chromosome <- c()
targetStart <- c()
targetEnd <- c()
ampliconName <- c()

for( amplicon_ID in unique(paa_pl_a_pl_p$amplicon_id))
{
  # target start
  tStart <- filter( paa_pl_a_pl_p, amplicon_id == amplicon_ID &  primer_strand == 'forward')$primer_end
  tStart <- tStart + 1
  targetStart <- c(targetStart, tStart)
  # target end
  tEnd <- filter( paa_pl_a_pl_p, amplicon_id == amplicon_ID &  primer_strand == 'reverse')$primer_start
  tEnd <- tEnd - 1
  targetEnd <- c(targetEnd, tEnd)
  
  # chromosome
  chr <- filter( paa_pl_a_pl_p, amplicon_id == amplicon_ID &  primer_strand == 'forward')$chromosome
  chromosome <- c(chromosome, chr)
  
  # amplicon name
  amplicon_name <- filter( paa_pl_a_pl_p, amplicon_id == amplicon_ID &  primer_strand == 'forward')$amplicon_name
  ampliconName <- c(ampliconName, amplicon_name)
}

targetTable <- data.frame( chromosome=chromosome, targetStart=targetStart, 
                           targetEnd=targetEnd, strand=rep(strand, length(chromosome)), ampliconName=ampliconName)

targetTable <- arrange(targetTable, chromosome, targetStart)

write.table(x=targetTable, file='targets_v1.0.txt',
            row.names = FALSE, col.names = FALSE, quote=FALSE, se='\t'
           )

