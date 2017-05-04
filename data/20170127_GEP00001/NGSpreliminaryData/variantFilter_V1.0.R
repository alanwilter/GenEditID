# Script for filtering indels based on various factors
# Author - Chandu
# date : 04-05-17

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

# Do we consider depth filter? #Ruben As long as the variant caller gives 'pass', I think we are OK.

# cut sites
# Need to modify script to get the cut sited from database
# currently hard coded in the script
# STAT3.2.in__STAT3.3.in__STAT3.4.in = 40497640
# STAT3.1.in = 40498700
# STAT3.3.off1= 125765667
tab <- read.delim( file='SLX-13775_variantsINDELS.csv', stringsAsFactors = F)  #Ruben. Add variantsSNV's

# clean col names
names(tab) <- gsub('\\.','_',names(tab) )
names(tab) <- gsub('__','_',names(tab) )
names(tab) <- sub('_$','', names(tab))

names(tab)[match( 'X5_context', names(tab) ) ] <- 'fivePrimeContext'
names(tab)[match( 'X3_context', names(tab) ) ] <- 'threePrimeContext'

# remove empty columns
tab <- select(tab, -Existing_variation : -GMAF, -PubMed) 

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
# Filter out AF > 90 and indel length > 3 #Ruben. Instead of filter out, weigh down. 
                                          # These could still be valid clones, although we don't have confidence to call them.

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
