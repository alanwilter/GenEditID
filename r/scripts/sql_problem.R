sql <- 
  'select
p.geid as Plate,
concat(w.row, w.column) as Well,
g.name as Guide,
slc.sequencing_sample_name as Sample,
slc.sequencing_barcode as Barcode,
vr.consequence as Consequence,
vr.gene_id as SymbolGene,
vr.cdna_effect as CDNAEffect,
vr.protein_effect as ProteinEffect,
vr.codons as Codons,
vr.chromosome as Chromosome,
vr.position as Position,
vr.sequence_ref as Ref,
vr.sequence_alt as Alt,
vr.allele_fraction as AlleleFraction,
vr.depth as Depth,
vr.quality as Quality,
vr.amplicon as Amplicon,
vr.gene as Gene,
vr.exon as Exon,
vr.offset_from_primer_end as PrimerEndOffset,
vr.indel_length as IndelLength,
vr.sequence_alleles as Alleles,
vr.variant_type as Type
from
well w
left join sequencing_library_content slc on slc.well_id = w.id
inner join variant_result vr on vr.sequencing_library_content_id = slc.id
inner join experiment_layout el on w.experiment_layout_id = el.id
left join plate p on p.experiment_layout_id = el.id
inner join well_content wc on w.well_content_id = wc.id
left join guide_well_content_association gwca on gwca.well_content_id = wc.id
inner join guide g on gwca.guide_id = g.id'
  
  NGSdata <- dbGetQuery(db$con, sql) %>%
    filter(allelefraction > 0.1)         # rows with allelefraction < 0.1 are likely PCR or sequencing errors
  
  # TO FIX!
  # Plate should be filtered to yield only "NGS plates"
  NGSdata <- filter(NGSdata, grepl('incu', plate)) # temporary workaround until we fix the plate_NGS problems
  
  NGSdata[!mapply(grepl, NGSdata$guide, NGSdata$amplicon),]