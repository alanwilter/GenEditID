# Planning - What's next? Our TODO list!

## layout data template
file `data/templatesYYYYMMDD_GEPXXXXX.xlsx`

#### Done
- ANNE (2) [x] we don't need 'barcode_size' in SequencingLibraryContent
  - remove it everywhere from model, template and loader
- ANNE (2) [x] in tab GuideMismatches: 'regions' would better be 'number_of_mismatches', and 'mismatches' should be 'number_of_offtargets'
- ANNE (2) [x] in Project: change 'institute' to 'affiliation'

## data model (model.py)
- ANNE (.) [ ] check that codes associated to classes are not in script but put back in class functions
- ANNE (?) [ ] Add donor table: sequence, start, end on forward and excision sequence

#### Done
- ANNE (2) [x] add 'cell_pool' to ExperimentLayout
- ANNE (2) [x] add (editable) comments on projects
- ANNE (2) [x] Add column in mutation_summary for found in both variant callers VH, V-, -H, or V?
- ANNE (2) [x] Add column in mutation_summary for frameshift yes/no
- ANNE (2) [x] add column in mutation_summary for score
- ANNE (2) [x] Add mismatches table for guides: is_coding_region, regions, mismatches
- ANNE (3) [x] Add 'project_type' ('knock-in', 'knock-out') to Project.
- ANNE (3) [x] update database diagram

## loader
- ANNE (3) [ ] Add checks and error messages to ensure layouts are correct (e.g. in ExperimentLayout, a .content_type "wild-type" cannot have a guide associated), and a .content_type "empty" should not have anything else associated other than well position, a "sample" should have all values (cell_line, clone, guide, replicate).
- [ ] In the loader, *ProteinAbundanceLoader* class, change file type from '.csv' to '.txt' (or it looks like extension is not used in the loader, only sep?). It's a tab-delimited file anyway, so this way we can keep it consistent with the growth data .txt extension. I am doing the documentation of protein and growth file formats and noting this down. The extension does not matter for the loader, it only needs to be the same format. We use the csv library to load the file and specify the delimiter to be tab. If we were renaming these files to .txt, we will need to update the loading scripts used to populate the database in `shell/` directory.
- [ ] In the loader, it may be better to calculate 'hour' from the timestamp. I'll talk to the incucyte techs about this, so don't change anything just yet.

#### Done
- ANNE (3) [x] calculate mutation_summary score and load
- ANNE (3) [x] add values for variant_caller_presence
- ANNE (2) [x] add values for has_frameshift
- ANNE (2) [x] load guide mismatches

## sequencing analysis pipeline (later)
- CHANDU (3) [ ] genome coordinates. The user is using hg19 coordinates for primers and guides. We need a script to translate these coordinates to hg18.
- [ ] pipeline to work with human and mouse genomes
- [ ] Not urgent. primer design automation and loading data in DB (primer blast)

## pipeline automation (later)
- [ ] install pipeline in production as well as the software dependencies and databases for annotations
- [ ] check reference genomes in place
- [ ] create a space on the cluster for genome editing projects or elsewhere
- [ ] propose solution on full automation

## plots / analysis
- RUBEN (1) [ ] 96-well plate scatter plot (Ruben - need the scores and the slopes)
- [x] Indelranges (Chandu)
- CHANDU (2) [ ] Type of mutation bar plot (%of samples submitted to NGS vs x =[wt, ins, del, SNV])
- CHANDU (2) [ ] Heatmap [protein, NGS[[frameshift-frame1, frameshift-frame2, noframeshift], [offtarget], [zygosity]], growth slope]
- RUBEN (2) [ ] Growth slopes (we need to calculate the slopes!) - done except with colors
- [ ] Combined plot NGS + protein + slopes, color-coded for frameshifts

- [ ] Not urgent. Current analysis is per single project > update to process multiple projects
- [ ] Not urgent. Current score calculations are based on allele number == 2. However in the PTEN project 3 alleles are knocked, and in future projects multiple positions could be edited simultaneously - need to add an 'allele_number' field to the database and update the code

## pyramid webapp
- **project description**
  - [ ] Core Genome Editing comments section: box to add comments about the project, that should be updatable and stored on the database along the project
  - [ ] project global overview: tab 'Project' from the layout excel file
  - [ ] project detailed overview (on demand, clicking somewhere): tabs 'Project', 'Target', 'Guide' and 'Guidemismatches' from the layout excel file
- **global data exploration**, plots for:
  - RICH [ ] Protein abundance
  - RICH [ ] Growth curves
  - [ ] Growth slopes
  - NGS exploratory plots
    - ANNE (1) [ ] % Zygosities (‘homozygous’, ‘heterozygous’…)
    - [ ] % Alleles (e.g. C/CATG, CTAA/C)
    - [ ] % Distances from cut site
    - [ ] % Indel lengths
    - [ ] % Type of mutation (‘wt’, ‘insertion’, ‘deletion’, ‘SNV’)
    - [ ] % Type of variant (‘wt’, ‘frameshift’, ‘inframe’...)
    - [ ] INDEL ranges (visual location of alleles)
  - [ ] Combined plots (NGS + protein + growth slopes)
  - [ ] Heatmap
- **sample details**
  - [ ] plot: 96-well plates for score-based clone selection
  - RICH [ ] Data table per project
  - [ ] Per-sample plot of INDEL ranges
  - LATER [ ] Visualisation of reads (.bam files) in IGV browser (external) (ideally a link to the .bam that opens an IGV server, or an IGV installed locally on the user's computer)
- **user report** for selected samples
  - need a way to select samples and get data table (and plots) only related to these
  - [ ] project global overview (tab 'Project' from the layout excel file)
  - [ ] plot INDEL ranges
  - [ ] Combined plot (clone vs controls)
  - Data table
    - [ ] per project and plate layout (i.e. GEP00001_01 and no GEP00001_01_incu)
    - [ ] output data as soon as we have input (e.g. if there is protein data, show it regardless of presence of growth and NGS data). Ideally we could click on the 96-well plate and have the clone highlighted in this table, but this is not priority
    - [ ] columns of the data table:
      - Plate
      - Well
      - Sample name
      - Fluidigm barcode
      - Score (output NaN if no NGS data available)
      - Protein abundance (800/700) relative to wild type control (this needs to be calculated. Output NaN if no protein data available)
      - Protein abundance (800/700) relative to knock-out control (this needs to be calculated)
      - (Growth slope relative to wild type control (this needs to be calculated. Output NaN if no growth data available))
      - (Growth slope relative to knock-out control (this needs to be calculated))
      (these below, columns from the variant output. Output NaN if no NGS data available )
      - Variant type/consequence
      - Symbol (Gene ID)
      - IGV link
      - Allele fraction
      - Alleles
  - [ ] Comments (this section is independent to the one described in 'project')
- [ ] generate final user report

## problems to solve
- [x] in plot_typeofvariants.py and plot_distances.py, weird result: for guide STAT3.1 I get a -91 indel_length with haplotypecaller, but according to the excel file that should be STAT3.3 instead (sample GE-P1B5-C).
- [x] is presence of offtargets considered in calculation of zygosities?
- [ ] revise queries, filtering by project might not be working correctly (we can use sample names as tracking system)
- [ ] very large indels detected (90 bp!)
- [ ] issue with *shiny app*: position_dodge messes up hover event data (Ruben is on it, raised issue at github plot.ly)
- CHANDU (1) [x] haplotypecaller is giving two merged values for the same sample onto one row separated by comma
- [ ] unable to fit a model on growth plots to extract curve (Discussion with Dominique on Thu 13 July)
