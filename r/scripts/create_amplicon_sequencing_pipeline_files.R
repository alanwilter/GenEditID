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
# target_start = forward_primer_end + 1
# target_end  = reverse_primer_start - 1

# --------------------------------------------------------------------------------
# Questions?
# --------------------------------------------------------------------------------
# (1) Where does amplicon seq pipeline use strand information?
# Not clear, needs to discuss with Matt
# Currently all the coordinates are forward strand based
# (2) Should target name be different from amplicon name?
# No, both amplicon and target names should be identical, therefore using amplicon names for both

# --------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------
# load data from database using amplicon query
loadData <- function(db_connection, slxid, projectid)
{
  sql <- paste0(
    "select distinct
        sequencing_library.slxid,
        cast(amplicon.chromosome as integer) as chr_number,
        -- amplicon
        'chr' || amplicon.chromosome as amplicon_chr,
        amplicon.start as amplicon_start,
        amplicon.end as amplicon_end,
        '+' as amplicon_strand,
        genome.assembly || '_chr' || amplicon.chromosome || '_' || amplicon.start as amplicon_name,
        -- target
        'chr' || amplicon.chromosome as target_chr,
        forward_primer.end + 1 as target_start,
        reverse_primer.start - 1 as target_end,
        '+' as target_strand,
        genome.assembly || '_chr' || amplicon.chromosome || '_' || forward_primer.end as target_name
    from sequencing_library,
        sequencing_library_content,
        well,
        well_content,
        guide_well_content_association,
        guide,
        target,
        project,
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
        and guide.target_id=target.id
        and target.project_id=project.id
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
        and sequencing_library.slxid='", slxid, 
      "'and project.geid='", projectid, "'
    order by chr_number, amplicon_start")

  query <- as.data.frame(dbGetQuery(db_connection, sql))
  query
}

# --------------------------------------------------------------------------------
# Main code
# --------------------------------------------------------------------------------
driver <- dbDriver("PostgreSQL")
conn <- dbConnect(driver, user="gene", password="gene", host="bioinf-ge001.cri.camres.org", port=5432, dbname="geneediting")
projectid <- 'GEP00002'
slxid <- 'SLX-13774'
data <- loadData(conn, slxid, projectid)
print(data)

# amplicon coordinates
write.table(x=select(data, contains('amplicon')), file='amplicons.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, se='\t')

# target coordinates
write.table(x=select(data, target_chr, target_start, target_end, target_strand, amplicon_name),
            file='targets.txt', row.names=FALSE, col.names=FALSE, quote=FALSE, se='\t')
