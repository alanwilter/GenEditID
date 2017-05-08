# Planning meetings

## Data in database

- store NGS results

## Variant calling methods

- three methods, two methods used so far and make union
- taking one method to start with and thinking of unions and weight later

## Plots: Results to show

### Static plots first
- plots only against a specific project
- plot growth curves
- plot protein against growth slopes - add smoothing to the calculation of the slope
- show the right clone by giving a score and displaying on a plate view (red/green)
  - allele frequency - filtering, threshold 0.15 and then take into account for calculating score
  - weight down the ones with off-target
  - real mutation
  - location of indel - if location is the same of the target (taking frameshift into account)

### Interactive plots next

## Plan

- Anne will load the NGS results into the database
- Rich will do the plate drawing based on the plate layout project from github https://github.com/crukci-bioinformatics/PlateLayout
- Ruben will do the smoothing of the slope [r/scripts/genome_editing.r](../r/scripts/genome_editing.r)
- Chandu will do the calculation of the NGS score

[Database installation/information](postgres.md)

# TODO list

- [] Need to modify variantFilter_V1.0.R to get the cut site from database
- [] load data from both indels and snvs into variantFilter_V1.0.R
- [] issue with shiny app: position_dodge messes up hover event data. Ruben is on it, raised issue at github plot.ly.
- [] Current score calculations are based on allele number == 2. However in the PTEN project 3 alleles are knocked, and in future projects multiple positions could be edited simultaneously. We need to add an 'allele_number' field to the database.
- [] There are large frameshifts missing in either haplotype or Vardict callers, so it could be that we need to merge results from both?
