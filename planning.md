# Planning meetings

## Data in database
- store NGS results

[Database installation/information](postgres.md)


## Variant calling methods

- three methods, two methods used so far and make union
- taking one method to start with and thinking of unions and weight later

## Plots: Results to show

### Static plots first
- plots only against a specific project
- plot growth curves
- plot protein against growth slopes - add smoothing to the calculation of the slope
- show the right clone by giving a score and displaying on a plate view (red/green) using the plate layout project from github https://github.com/crukci-bioinformatics/PlateLayout
  - allele frequency - filtering, threshold 0.15 and then take into account for calculating score
  - weight down the ones with off-target
  - real mutation
  - location of indel - if location is the same of the target (taking frameshift into account)

### Interactive plots next

## Today's plan
- [ ] Rich: load data from db in shiny app
- [ ] Ruben: put Rich plot in shiny app
- [ ] Chandu: modify script to read from db
- [ ] Chandu: move script
- [x] Anne: git move?
- [x] Anne: clean python loader code
- [ ] Anne: add mismatch table
- [ ] Anne: look at python webapp options and framework


## What's next? Our TODO list...
- [ ] need to modify variantFilter_V1.0.R to get the cut site from database
- [ ] load data from both indels and snvs into variantFilter_V1.0.R
- [ ] issue with shiny app: position_dodge messes up hover event data (Ruben is on it, raised issue at github plot.ly)
- [ ] Current score calculations are based on allele number == 2. However in the PTEN project 3 alleles are knocked, and in future projects multiple positions could be edited simultaneously - need to add an 'allele_number' field to the database
- [ ] There are large frameshifts missing in either haplotype or Vardict callers, so it could be that we need to merge results from both?

**Following Friday 5th May 2017 3 month review:**
- primer design automation and loading data in DB (primer blast)
- guide vs variants score given by the website (Ruben has code)
- base content of variants and coordinates
- analyse different projects together
- web app
