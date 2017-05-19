library( 'dplyr')
library( 'RPostgreSQL')

# --------------------------------------------------------------------------------
# Definitions
# --------------------------------------------------------------------------------
# <--------- Amplicon ---------->
# ===============================
# forward primer   reverse primer
# >>>>>>>                 <<<<<<<
#        <---- Target --->

# --------------------------------------------------------------------------------
# Questions?
# --------------------------------------------------------------------------------
# (1) Where does amplicon seq pipeline use strand information?
# Not clear, needs to discuss with Matt
# Currently all the coordinates are forward strand based
# (2) Should target name be different from amplicon name?

# --------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------
# load data from database using amplicon query
loadData <- function(db_connection, slxid) 
{
  sql <- paste0(
    "select distinct 
        sequencing_library.slxid,
        -- amplicon
        'chr' || amplicon.chromosome as amplicon_chr, 
        amplicon.start as amplicon_start, 
        amplicon.end as amplicon_end,
        case when amplicon_selection.guide_strand = 'forward' then '+' else '-' end as amplicon_strand,
        genome.assembly || '_chr' || amplicon.chromosome || '_' || amplicon.start as amplicon_name, 
        -- target
        'chr' || amplicon.chromosome as target_chr, 
        forward_primer.end as target_start, 
        reverse_primer.start as target_end,
        case when amplicon_selection.guide_strand = 'forward' then '+' else '-' end as target_strand,
        genome.assembly || '_chr' || amplicon.chromosome || '_' || forward_primer.end as target_name
    from sequencing_library, 
        sequencing_library_content, 
        well, 
        well_content, 
        guide_well_content_association, 
        guide, 
        amplicon_selection, 
        amplicon, 
        genome, 
        primer_amplicon_association as forward_association, 
        primer_amplicon_association as reverse_association, 
        primer as forward_primer, 
        primer as reverse_primer
    where sequencing_library_content.sequencing_library_id=sequencing_library.id
        and sequencing_library_content.well_id=well.id
        and well.well_content_id=well_content.id
        and guide_well_content_association.guide_id=guide.id
        and guide_well_content_association.well_content_id=well_content.id
        and amplicon_selection.guide_id=guide.id
        and amplicon_selection.amplicon_id=amplicon.id
        and amplicon.genome_id=genome.id
        and forward_association.amplicon_id=amplicon.id
        and forward_association.primer_id=forward_primer.id
        and forward_primer.strand='forward'
        and reverse_association.amplicon_id=amplicon.id
        and reverse_association.primer_id=reverse_primer.id
        and reverse_primer.strand='reverse'
        and sequencing_library.slxid='", slxid, "'")

  query <- as.data.frame(dbGetQuery(db_connection, sql))
  query
}

# --------------------------------------------------------------------------------
# Main code
# --------------------------------------------------------------------------------
driver <- dbDriver("PostgreSQL")
conn <- dbConnect(driver, user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")
slxid <- 'SLX-13775'
data <- loadData(conn, slxid)
print(data)

# amplicon coordinates
amplicon_coordinates <- select(data, contains( 'amplicon') ) %>% 
  arrange( amplicon_chr, amplicon_start)
write.table(x=amplicon_coordinates, file='amplicons.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, se='\t')

# target coordinates
# target_start = forward_primer_end + 1
# target_end  = reverse_primer_start - 1
target_coordinates <- select(data, contains( 'target') ) %>% 
  mutate( target_start = target_start + 1, target_end=target_end-1) %>% 
  arrange( target_chr, target_start)
write.table(x=target_coordinates, file='targets.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, se='\t')

