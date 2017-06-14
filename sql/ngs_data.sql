with guide_list as (
   select
      guide_well_content_association.well_content_id,
      string_agg(guide.name, ', ') as guide_names
   from
      guide_well_content_association
      inner join guide on guide_well_content_association.guide_id = guide.id
   group by guide_well_content_association.well_content_id
)

select
   experiment_layout.geid as experiment_layout_geid,
   concat(well.row, well.column) as well,
   cell_line.name as cell_line_name,
   clone.name as clone_name,
   guide_list.guide_names,
   well_content.is_control,
   well_content.content_type,
   sequencing_library_content.dna_source,
   sequencing_library_content.sequencing_sample_name,
   sequencing_library_content.sequencing_barcode,
   sequencing_library.slxid,
   variant_result.consequence,
   variant_result.gene_id,
   variant_result.cdna_effect,
   variant_result.protein_effect,
   variant_result.codons,
   variant_result.chromosome,
   variant_result.position,
   variant_result.sequence_ref,
   variant_result.sequence_alt,
   variant_result.allele_fraction,
   variant_result.depth,
   variant_result.quality,
   variant_result.amplicon,
   variant_result.gene,
   variant_result.exon,
   variant_result.offset_from_primer_end,
   variant_result.indel_length,
   variant_result.sequence_alleles,
   variant_result.variant_type,
   variant_result.ge_score
from
   well
   inner join experiment_layout on well.experiment_layout_id = experiment_layout.id
   left join well_content on well.well_content_id = well_content.id
   left join guide_list on guide_list.well_content_id = well_content.id
   left join clone on well_content.clone_id = clone.id
   left join cell_line on clone.cell_line_id = cell_line.id
   left join sequencing_library_content on sequencing_library_content.well_id = well.id
   inner join sequencing_library on sequencing_library_content.sequencing_library_id = sequencing_library.id
   inner join variant_result on variant_result.sequencing_library_content_id = sequencing_library_content.id
where
   variant_result.allele_fraction > 0.1
