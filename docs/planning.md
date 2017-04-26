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
