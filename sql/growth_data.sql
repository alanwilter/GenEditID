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
   growth.hours,
   growth.confluence_percentage
from
   well
   inner join experiment_layout on well.experiment_layout_id = experiment_layout.id
   left join well_content on well.well_content_id = well_content.id
   left join guide_list on guide_list.well_content_id = well_content.id
   left join clone on well_content.clone_id = clone.id
   left join cell_line on clone.cell_line_id = cell_line.id
   left join growth on growth.well_id = well.id
   left join sequencing_library_content on sequencing_library_content.well_id = well.id
