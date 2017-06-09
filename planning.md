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
- [x] Rich: load data from db in shiny app
- [x] Ruben: put Rich plot in shiny app
- [x] Chandu: modify script to read from db
- [x] Chandu: move script
- [x] Anne: git move?
- [x] Anne: clean python loader code
- [ ] Anne: add mismatch table
- [x] Anne: look at python webapp options and framework


## What's next? Our TODO list...
- [ ] 1. tidy up scripts
  - [ ] 1b.layout_plot.R to use data from all plates
- [ ] 2. meet with bioformatics teams from Sanger and AZ
- [ ] 3. current analysis is per single project > update to process multiple projects
- [ ] 4. prepare a knock-in donor table to include in database
- [ ] 5. pipeline to work with human and mouse genomes
- [x] 6. need to modify variantFilter_V1.0.R to get the cut site from database
- [x] 7. load data from both indels and snvs into variantFilter_V1.0.R
- [ ] 8. generate final user report from app (html report?)
- [ ] 9. issue with shiny app: position_dodge messes up hover event data (Ruben is on it, raised issue at github plot.ly)
- [ ] 10. Current score calculations are based on allele number == 2. However in the PTEN project 3 alleles are knocked, and in future projects multiple positions could be edited simultaneously - need to add an 'allele_number' field to the database
- [ ] 11. There are large frameshifts missing in either haplotype or Vardict callers, so it could be that we need to merge results from both?
- [ ] 12. currently, one sql query per dataset (NGS, growth, protein). Perhaps one query for all (that allows for missing datasets)

**Following Friday 5th May 2017 3 month review:**
- primer design automation and loading data in DB (primer blast)
- guide vs variants score given by the website (Ruben has code)
- base content of variants and coordinates
- analyse different projects together
- web app

## Problems to solve
- [ ] P1. In NGSplot.R, fun.NGS_readDB() gives a table with NGS data. In field 'plate', we get plates of all types (NGS, PCR, dilution...) and consequently duplicated barcodes. I have pushed a corrected version of the spreadsheet tab 'Plate' (20170127_GEP00001_problemPlate.csv, *not* updated in the Excel file), where all layouts are connected to their respective plate barcodes, so NGSdata can be filtered to yield only NGS plates.
- [ ] P2. from tab Guide add 'target_name' to fun.NGS_readDB() in NGSplot.R. This is needed to calculate zygosities in target and off-targets independently within fun.NGS_exploratory()
- [ ] P3. from tab AmpliconSelection add 'guide_location' to fun.NGS_readDB() in NGSplot.R. This is needed in NGSplot.R, line 192. In geom_vline(xintercept = guide_location, color = 'steelblue', linetype = "dotted"), 'xintercept' should be guide cut site (guide_location), dependent on guide and supplied from the database. Currently only valid for project1.
- [ ] P4. problem with NGS sql query (or database itself?), with incorrect assignment guide/amplicon. See result of running sql_problem.R (e.g. first line, guide 3.1 has amplicon  STAT3.2.in__STAT3.3.in__STAT3.4.in)
- [ ] P5. in growth_plot.R, fun.growth_readDB(db) gives a tbl with NA's in target_id and guide_id ('guide' information is needed to sort the plots)
- [x] P6. in growth_plot.R, something happened to the data, for some reason now we can't calculate slopes! (no errors, but slopes don't work). This could be related to the incorrect create_classifier() output (see that e.g. "MCF7 clone3 STAT3.3 normalisation" should not exist). This odd output applies also to the protein data query. In both growth_plot.R and abundance_plot.R we are using the dplyr queries. (@pajanne I do believe it is now sorted after fixing a bug in the loader)
  ```
    unique(clone_growth_data$Content)
     [1] "MCF7 clone3 normalisation"         "MCF7 clone3 STAT3.3"               "MCF7 clone3 STAT3.3 normalisation"
     [4] "MCF7 clone3 STAT3.4"               "MCF7 clone3 STAT3.4 normalisation" "MCF7 clone3 STAT3.1"              
     [7] "MCF7 clone3 STAT3.2"               "MCF7 clone3 STAT3.3 background"    "MCF7 clone3 STAT3.3 knock-out"    
    [10] "MCF7 clone3 STAT3.3 wild-type"
  ```  
- [ ] P7. in growth.data <- fun.growth_readDB(db), column 'barcode' is incorrect and shows the same value as column 'geid'
- [ ] P8. regarding point 12 in TODO list, currently we are merging growth, protein and NGS to generate a combined plot. Merging by plate and well is not possible because they are not unique (we have cells or gDNA coming from the same well), so the ideal merger would be 'sample', but it's not present in all datasets. Instead of trying to fix this, would it be better to go for a single query, where the query gives a value for presence/absence of a certain dataset (growth, protein, NGS) that tells us which plots we can make?
