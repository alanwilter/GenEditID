tab <- read.delim( file='SLX-13775_variantsSNVs.csv', stringsAsFactors = F)

# clean col names
names(tab) <- gsub('\\.','_',names(tab) )
names(tab) <- gsub('__','_',names(tab) )
names(tab) <- sub('_$','', names(tab))

names(tab)[match( 'X5_context', names(tab) ) ] <- 'fivePrimeContext'
names(tab)[match( 'X3_context', names(tab) ) ] <- 'threePrimeContext'

# remove empty columns
tab <- select(tab, -Existing_variation : -GMAF, -PubMed) 
tab <- mutate(tab, OnTarget=rep(TRUE, nrow(tab)),
              cutSite=rep(40497640, nrow(tab)),
              distance=rep(0, nrow(tab)),
              disScore=rep(1, nrow(tab))
              )

write.csv( x=tab, file='SLX-13775_variantsSNVs_filtered.csv', row.names=F)
