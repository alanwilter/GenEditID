# Script for filtering indels based on various factors
# Author - Chandu
# date : 04-05-17

library( dplyr)
library( ggplot2)


cleanColNames <- function( tab, varType, varCaller ){
  # clean col names
  names(tab) <- gsub('\\.','_',names(tab) )
  names(tab) <- gsub('__','_',names(tab) )
  names(tab) <- sub('_$','', names(tab))
  
  names(tab)[match( 'X5_context', names(tab) ) ] <- 'fivePrimeContext'
  names(tab)[match( 'X3_context', names(tab) ) ] <- 'threePrimeContext'
  tab <- mutate( tab, varType=rep(varType, nrow(tab)),
                 varCaller=rep(varCaller, nrow(tab))
                 )
  return(tab)
  
}

# Need to add additional column indicating on or off targets
# filter based on alle frequency
# Filter all the sites <0.15 ( equvivalent of 6 alleles 100/6 ~ .15%)
# Filter all the alleles show high frequency ( > 90),
# which indicates there are two alleles 
# 0.4 -0.6 allele frequency indicates heterozygous

# It may depend on the length of the indel, longer the indel, rarer the homozygous probability
# need to filter based on the length, but what is the cut off?

# Filter based on the effect

# Filter based the relative location of the mutation from ther target site

# Where to get the exact target position?

# Do we consider depth filter? #Ruben As long as the variant caller gives 'pass', I think we are OK.

# cut sites
# Need to modify script to get the cut sited from database
# currently hard coded in the script
# STAT3.2.in__STAT3.3.in__STAT3.4.in = 40497640
# STAT3.1.in = 40498700
# STAT3.3.off1= 125765667

# read cut site info from DB

# Connect to DB
db <- src_postgres(user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting" )

# Data is spread in two tables
# amplicon_selection
# amplicon

# cut site = amplicon_selection.guide_location
# Amplicon Name = amplicon.name

# read amplicon_selection table
amplicon_selection <- tbl( db, 'amplicon_selection') %>% 
  select( amplicon_id, guide_location) %>% 
  as.data.frame()


# read amplicon table
amplicon <- tbl( db, 'amplicon' ) %>% 
  select( id, name) %>% 
  rename( amplicon_id=id) %>% 
  as.data.frame()

# merge two tables
cutInfoTab <- left_join(amplicon_selection, amplicon, by='amplicon_id' )

tab <- read.delim( file='../../data/20170127_GEP00001/NGSpreliminaryData/SLX-13775_variantsINDELS.csv', stringsAsFactors = F)  #Ruben. Add variantsSNV's

# clean column names
tab <- cleanColNames( tab,varType='indel', varCaller='VarDict')

# add on/off target info
tab <- mutate(tab, OnTarget=rep('-', nrow(tab)))
tab$OnTarget[ grep( 'off', tab$Amplicon)] <- 'FALSE'
tab <- mutate( tab,OnTarget=ifelse(OnTarget == 'FALSE', FALSE, TRUE) )


# add cut site info
tab <- mutate(tab, cutSite=ifelse(Amplicon =='STAT3.2.in__STAT3.3.in__STAT3.4.in', 40497640,  
                                  ifelse( Amplicon == 'STAT3.1.in', 40498700, 125765667)))
tab <- mutate(tab, distance=abs(Position - cutSite))


# Distance   Score
# 0          4  # on target
# 1-5        3
# 6-10       2
# > 10       1

# disScore = distance based score

tab <- mutate( tab, disScore=rep( 0, nrow(tab))) 
tab[ tab$distance ==0, 'disScore'] <- 4
tab[ (tab$distance >=1 & tab$distance < 6), 'disScore']  <- 3
tab[ (tab$distance >=6 & tab$distance < 11), 'disScore']  <- 2
tab[ tab$distance > 10, 'disScore']  <- 1

# Wt /Het/ Home
# Classification based on allele frequency
# Het = 40-60%

# Filter out AF < 0.15 ( probably wt or false positives)
# Filter out AF > 90 and indel length > 3 #Ruben. Instead of filter out, weigh down. 
                                          # These could still be valid clones, although we don't have confidence to call them.

# Allele frequency distribution plot
p <- ggplot( data=tab, aes( x=Allele_fraction)) + 
       geom_histogram( bins = 11, colour='orange', fill='white') + 
       ggtitle( 'Allele fraction distribution') +
       xlab( 'Allele Fraction')
ggsave('totalAlleleFractionDistribution.pdf', p)


# Allele fraction distribution > 0.9
p <- ggplot( data=filter(tab, Allele_fraction > 0.9), aes( x=Allele_fraction)) +
  geom_histogram( bins = 11, colour='orange', fill='white') +
  ggtitle( 'Allele fraction distribution > 0.9') +
  xlab( 'Allele Fraction')
ggsave( 'higherEndAlleleFractionDistribution.pdf', p)

# Allele fraction distribution < 0.15

p <- ggplot( data=filter(tab, Allele_fraction < 0.15), aes( x=Allele_fraction)) +
  geom_histogram( bins = 11, colour='orange', fill='white') +
  ggtitle( 'Allele fraction distribution < 0.15') +
  xlab( 'Allele Fraction')
ggsave( 'lowerEndAlleleFractionDistribution.pdf', p)


# Total indel length distribution 
p <- ggplot( data=tab, aes( x=Indel_length)) +
  geom_histogram( colour='orange', fill='white') +
  ggtitle( 'Indel Length Distribution') +
  xlab( 'Indel Length')
ggsave( 'totalIndelLengthDistribution.pdf', p)

#  indel length distribution > 0.9 allele fraction
p <- ggplot( data=filter(tab, Allele_fraction > 0.90 ), aes( x=Indel_length)) +
  geom_histogram( colour='orange', fill='white') +
  ggtitle( 'Indel Length Distribution > 0.9 Allele Fraction') +
  xlab( 'Indel Length')
ggsave( 'indelLengthDistributionAtHigherAlleFraction.pdf', p)



  
  
tab <- filter(tab, Allele_fraction > 0.15)

tab <- filter(tab, !(Allele_fraction > 0.9 & abs(Indel_length) > 3  ))

# filter based on consequence #Ruben. Weigh down instead. Inframe_deletion could generate a phenotype if it is on a key position 
                              #in the protein.
                              #Very careful with filtering, because values in 'sample' are not unique. If you filter out 'intron,non_coding_transcript'
                              #you will think that GE-P1B4-C is a good HET line because of alleleFrequency and OnTarget == TRUE. However it should be heavily weighed down because
                              #it's showing an off-target. When filtering, do it with effect on 'sample' groups, not on rows.

tab <- filter(tab, Variant_type_consequence == 'frameshift' | 
                Variant_type_consequence == 'stop_gained,frameshift')

# Filter based on on/off-target
tab <- filter(tab, OnTarget == TRUE)

write.csv( x=tab, file='SLX-13775_variantsINDELS_filtered.csv', row.names=F)


#Ruben. Suggested scoring system (the higher the better): what do you think, does it look OK from your experience? If not, feel free to modify
  #so we can test it later
 
  #Variant.type.consequence:
   # "frameshift": 5
   # "inframe_deletion": 1
   # "intron,non_coding_transcript": 0 #note that this will super-rarely happen for targets, because we would not normally target an intron.
                                       #so this almost always going to be the case of an off-target.
   # "stop_gained,frameshift": 6
 
 #Allele.fraction
   # Filter out A.f < 0.15
   # A.f > 0.4 & A.f < 0.6: score 5 (however we are assuming that we have 2 alleles here, and in project 2 we have three alleles.
                                    #either we are more lax, like >0.3 & <0.7, or we input the number of alleles from somewhere else)
   # A.f <=0.4 & A.f >=0.6 & A.f. <=0.9: score 1 (it could be real, but allele proportions are odd. Cross-contamination?)
   # A.f >0.9 & distance <4: score 4 (it could be real)
   # A.f >0.9 & distance >=4: score 2 (it could be real, but not enough data to know)

  #OnTarget (per sample group)
   # absence/presence of off-target in sample group: score 10/0


  #distance
   # d < 4: score 5
   # d >= 4: score 2

  #from protein dataset: use these scores only for KO projects (for knock-ins it would not make sense)
   # ratio protein clone/KO < 0.2: score 5 (probably KO)
   # ratio >=0.2 & <= 0.7: score 3 (probably HET wt/mutant)
   # ratio >0.7: score 0 (protein too high, even in the case of three alleles with one OK and the others WT)

  #from growth slope dataset:
   # Not needed, we can just sort the results table based on growth rate.
