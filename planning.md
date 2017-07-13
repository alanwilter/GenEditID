# Planning - What's next? Our TODO list!

## layout data template
file `data/templatesYYYYMMDD_GEPXXXXX.xlsx`

- [ ] we don't need 'barcode_size' in SequencingLibraryContent
- [ ] in tab GuideMismatches: 'regions' would better be 'number_of_mismatches', and 'mismatches' should be 'number_of_offtargets'
- [ ] in Project: change 'institute' to 'affiliation'

## data model (model.py)
- [ ] add 'pool' to ExperimentLayout
- [ ] prepare a knock-in donor table to include in database (not urgent)

## sequencing analysis pipeline
- [ ] genome coordinates. The user is using hg19 coordinates for primers and guides. We need a script to translate these coordinates to hg18.
- [ ] pipeline to work with human and mouse genomes
- [ ] Not urgent. primer design automation and loading data in DB (primer blast)

## plots / analysis
- [ ] 96-well plate scatter plot (Ruben)
- [ ] Indelranges (Chandu)
- [ ] Type of mutation bar plot (%of samples submitted to NGS vs x =[wt, ins, del, SNV])
- [ ] Heatmap [protein, NGS[[frameshift-frame1, frameshift-frame2, noframeshift], [offtarget], [zygosity]], growth slope]
- [ ] Growth slopes (we need to calculate the slopes!)
- [ ] Combined plot NGS + protein + slopes, color-coded for frameshifts

- [ ] Not urgent. Current analysis is per single project > update to process multiple projects
- [ ] Not urgent. Current score calculations are based on allele number == 2. However in the PTEN project 3 alleles are knocked, and in future projects multiple positions could be edited simultaneously - need to add an 'allele_number' field to the database and update the code
- [ ] Not urgent. currently, one sql query per dataset (NGS, growth, protein). Perhaps one query for all (that allows for missing datasets)

## pyramid webapp
- **project description**
  - [ ] Core Genome Editing comments section: box to add comments about the project, that should be updatable and stored on the database along the project
  - [ ] project global overview: tab 'Project' from the layout excel file
  - [ ] project detailed overview (on demand, clicking somewhere): tabs 'Project', 'Target', 'Guide' and 'Guidemismatches' from the layout excel file
- **global data exploration**, plots for:
  - [ ] Protein abundance
  - [ ] Growth curves
  - [ ] Growth slopes
  - NGS exploratory plots
    - [ ] % Zygosities (‘homozygous’, ‘heterozygous’…)
    - [ ] % Alleles (e.g. C/CATG, CTAA/C)
    - [ ] % Distances from cut site
    - [ ] % Type of mutation (‘wt’, ‘insertion’, ‘deletion’, ‘SNV’)
    - [ ] % Type of variant (‘wt’, ‘frameshift’, ‘inframe’...)
    - [ ] INDEL ranges (visual location of alleles)
  - [ ] Combined plots (NGS + protein + growth slopes)
- **clone details**
  - [ ] plot: 96-well plates for score-based clone selection
  - [ ] Data table
  - [ ] Per-clone plot of INDEL ranges
  - [ ] Visualisation of reads (.bam files) in IGV browser (external) (ideally a link to the .bam that opens an IGV server, or an IGV installed locally on the user's computer)
- **user report** for selected clones
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
      - Confidence
      - Allele fraction
      - Alleles
  - [ ] Comments (this section is independent to the one described in 'project')
- [ ] generate final user report

## problems to solve
- [x] in plot_typeofvariants.py and plot_distances.py, weird result: for guide STAT3.1 I get a -91 indel_length with haplotypecaller, but according to the excel file that should be STAT3.3 instead (sample GE-P1B5-C).
- [ ] is presence of offtargets considered in calculation of zygosities?
- [ ] very large indels detected (90 bp!)
- [ ] issue with *shiny app*: position_dodge messes up hover event data (Ruben is on it, raised issue at github plot.ly)
- [ ] haplotypecaller is giving two merged values for the same sample onto one row separated by comma
- [ ] unable to fit a model on growth plots to extract curve (Discussion with Dominique on Thu 13 July)
