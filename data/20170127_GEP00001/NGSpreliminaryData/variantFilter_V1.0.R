# Script for filtering indels based on various factors
# Author - Chandu
# date : 02-05-17

library( dplyr)

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

# Do we consider depth filter?

# cut sites
# Need to modify script to get the cut sited from database
# currently hard coded in the script
# STAT3.2.in__STAT3.3.in__STAT3.4.in = 40497640
# STAT3.1.in = 40498700
# STAT3.3.off1= 125765667
tab <- read.csv( file='SLX-13775_variantsINDELS.csv', stringsAsFactors = F)
tab <- mutate(tab, OnTarget=rep('-', nrow(tab)))
tab$OnTarget[ grep( 'off', tab$Amplicon)] <- 'FALSE'
tab <- mutate( tab,OnTarget=ifelse(OnTarget == 'FALSE', FALSE, TRUE) )

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
# Filter out AF > 90 and indel length > 3

tab <- filter(tab, Allele.fraction > 0.15) 

tab <- filter(tab, !(Allele.fraction > 0.9 & abs(Indel.length) > 3  ))

# filter based on consequence
tab <- filter(tab, Variant.type.consequence == 'frameshift' | 
                Variant.type.consequence == 'stop_gained,frameshift')

# Filter based on on/off-target
tab <- filter(tab, OnTarget == TRUE)

write.csv( x=tab, file='SLX-13775_variantsINDELS_filtered.csv', row.names=F)



