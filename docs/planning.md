# Genome Editing Planning - What's next? Our TODO list!

## layout data template
file `data/templates/GEPXXXXX.xlsx` with documentation explaining each sheets and each columns.


## data model (model.py)
- [ ] add donor table: sequence, start, end on forward and excision sequence
- [ ] add documentation to model.py (e.g. we have a relationship one to many in well to abundances. This is because you can take repeated measurements from the same well)
- [ ] revise relationships model.py (e.g. we have a relationship one to many in sequencing_library_content to mutation_summaries. In this case we can't obtain a second mutation_summaries for the same library, so it's more meaningful if relationship is 1:1).


## loader
- [ ] Add checks and error messages to ensure layouts are correct (e.g. in ExperimentLayout, a .content_type "wild-type" cannot have a guide associated), and a .content_type "empty" should not have anything else associated other than well position, a "sample" should have all values (cell_line, clone, guide, replicate).
- [ ] In the loader, it may be better to calculate 'hour' from the timestamp. I'll talk to the incucyte techs about this, so don't change anything just yet.
- [ ] calculate mutation_summary score and load. Currently score uses has_offtargets, zygosity and consequence. Add protein when the experiment is KO.
- [ ] The score uses has_frameshift. However, has_frameshift is True if at least one of the alleles has a frameshift. Make it True only when both alleles have a frameshift (the reason: if one is frameshift and the other inframe, the result is presence of protein and partial KO, but with the current score it's highly rated)


## sequencing analysis pipeline
- [ ] pipeline to work with human and mouse genomes
- [ ] primer design automation and loading data in DB (primer blast)


## automation
- [x] install pipeline in production as well as the software dependencies and databases for annotations
- [x] check reference genomes in place
- [ ] create a space on the cluster for genome editing projects or elsewhere
- [ ] propose solution on full automation


## plots / analysis
- [ ] Heatmap [protein, NGS[[frameshift-frame1, frameshift-frame2, noframeshift], [offtarget], [zygosity]], growth slope]
- [ ] Growth slopes (we need to calculate the slopes!) - done except with colors
- [ ] Combined plot NGS + protein + slopes, color-coded for frameshifts
- [ ] Added code to get indel structures from single variables. It is in ngsplotter.indelstructure_plot and needs to be put in the right place, connected to the table
- [ ] Not urgent. Current analysis is per single project > update to process multiple projects
- [ ] Not urgent. Current score calculations are based on allele number == 2. However in the PTEN project 3 alleles are knocked, and in future projects multiple positions could be edited simultaneously - need to add an 'allele_number' field to the database and update the code


## web app
- **global data exploration**, plots for:
  - [ ] Growth slopes
  - NGS exploratory plots
    - [ ] % Distances from cut site (to import)
    - [ ] % Type of mutation (‘wt’, ‘insertion’, ‘deletion’, ‘SNV’)
    - [ ] INDEL ranges (visual location of alleles)
    - [ ] ngsplotter.indelstructure_plot to add to plotter
  - [ ] Combined plots (NGS + protein + growth slopes)
  - [ ] Heatmap
- **sample details**
  - [ ] Per-sample plot of INDEL ranges
  - [ ] Visualisation of reads (.bam files) in IGV browser (external) (ideally a link to the .bam that opens an IGV server, or an IGV installed locally on the user's computer)
- **user report** for selected samples
  - need a way to select samples and get data table (and plots) only related to these
  - [ ] plot INDEL ranges
  - [ ] Combined plot (clone vs controls)
  - Data table
    - [ ] Ideally we could click on the 96-well plate and have the clone highlighted in this table, but this is not priority
    - columns of the data table:
      - [ ] Protein abundance (800/700) relative to wild type control (this needs to be calculated. Output NaN if no protein data available)
      - [ ] Protein abundance (800/700) relative to knock-out control (this needs to be calculated)
      - [ ] (Growth slope relative to wild type control (this needs to be calculated. Output NaN if no growth data available))
      - [ ] (Growth slope relative to knock-out control (this needs to be calculated))
      (these below, columns from the variant output. Output NaN if no NGS data available )
      - [ ] IGV link
  - [ ] Results Comments (this section is independent to the one described in 'project')
- [ ] generate final user report


## Problems to solve
- [ ] revise queries, filtering by project might not be working correctly (we can use sample names as tracking system)
- [ ] very large indels detected (90 bp!)
- [ ] issue with *shiny app*: position_dodge messes up hover event data (raised issue at github plot.ly)
- [x] haplotypecaller is giving two merged values for the same sample onto one row separated by comma
- [ ] unable to fit a model on growth plots to extract curve (Discussion with Dominique on Thu 13 July)
- [ ] make plots and data tables considering that we can have more than one dna_source. As it is currently, e.g. for project GEP00001, well.sequencing_library_contents[0].dna_source is 'fixed cells' and well.sequencing_library_contents[1].dna_source is 'gDNA'. We are selecting only fixed cells for simplification, but we need to show both (the reason to have gDNA and cells was to be able to compare the sequencing results from both). Also, in project one there is a sample (GE-P6B4-G) that is gDNA only (it was sent for sequencing only as gDNA, without a 'fixed cells' counterpart), no fixed cells, so when you do well.sequencing_library_contents[0].dna_source, it results in 'gDNA'!
- [ ] Revise protein data values in plot_96_wells, how the range looks like (we should see edges )


- **What's next?**
  - (1) Bioinformatics pipeline:
    - [x] new version of the amplicon sequencing pipeline after fixes on merged data per line
  - (2) Coordinate translation to switch between reference genomes
    - [x] use latest versions of genomes for everything e.g. Homo sapiens [GRCh38] or Mus musculus [GRCm38]
  - (3) Primer design: make it automatic - primer pair if failed nested primers
  - Use replicates - update calculation instead of using average
  - Question? Several guides in single well: how do you put in the layout file and load it?
  - Access control to only genome editing core
    - report with pdf or link to website
    - bam files
  - Load NGS results for knock-in project 3
  - Export selection of results not all results
  - Add growth slope calculation


- **What needs to be finished?**
  - Current web app
    - [x] Results table need to be updated with the right columns and one line per samples
    - [ ] Add result comments on project and add it on edit project
    - [ ] Finish scoring with protein ratio 800/700 vs KO average
    - [ ] Scale scoring system to 100%
    - [ ] Merge growth and abundance plots into one
    - [x] Check project type: knock-in and knock-out - it is in!
    - [x] Sample tab: update well name on plot to match the one in table
  - [ ] Update doc on loading data via webapp and ngs results manually
  - [ ] Document the script methods - 3 to 4 lines of comments in plotter files


- **How it went? What could we do better?**
  - R/Shiny to Python/Pyramid: is it normal process? nope but lots of learning
  - Agile way of working is great
  - Python, Plotly, database learning
  - Not spend that much time on the shiny
  - 2 times three months deadline has a negative impact on the project flow
  - Allocated working time to work together like two days per week
  - Great to have meeting every morning at 10am, very useful


## Catchup meeting - Wednesday 10th January 2018
**Main aim:** Project specific autopsies and moving forward

- Project specific autopsies
  - Plots with missing data on web app
  - Results table with missing columns on web app
- Moving forward
  - improve submitted form — check if columns could be locked
  - create standard location for BAM files
  - add BAM files in browser then results linked to hyperlink in IGV to take you to the right file
  - review scoring system: add threshold on number of reads
  - automatic submission to genomics
  - fully automated process: from submission to genomics to amplicon pipeline to upload to genome editing system
  - knock-in issue if insertion of 200 fragment — too long for alignment of 300 reads
    - build specific reference genome for each project when we have long knockin

- Next
  - link on the website to results excel table and also to the BAM files
  - column in the excel table so that the barcode number/clone number/well position can all be linked as they are in the submission form
  - summary of the wild type reads in this table too, for example FLD0001/A1 should be a parental wild type control, it would be good to see this to check for contamination in wells
  - if wild type reads are found within clones it would be good for this also to be listed so that we know that we have a heterozygous mutation
  - next project is a knock in project to run through the pipeline
    - over a 21bp region knock in 7 point mutations on both alleles
    - combining barcodes: screen the off targets for ATG5 in the same library


## End of project meeting - Thursday 5th Oct 2017
**Main aim:** Paper to publish: dry and wet - portable as much as possible: maybe docker

- Introduction
  - Is it useful what we currently have?
  - What is incomplete?
  - Smoothing the process at the NGS end
  - New projects to go through the system
  - Data presented in a more simplified way
  - Project3: knock-in build-in genome - why data not in?
  - Open system
  - Data analysis is a manual step
  - Amy point of contact from now on
  - NGS submission - library prep done
    - Next step: Submitting data into Genomics and within the system
    - Human genome: all knock-out

- Discussion
  - Data loader: [support Amy for filling spreadsheet]
    - support formula to load value?
    - amplicon and primer coordinates: checking tools?
      - guide RNA design - DeskGen
      - primer design
  - Primer design:
    - automatisation of primer blast: no downloadable version
      - currently primer3 but not the same as primer blast
      - maybe checking primer3 vs primer blast?
      - blat to check if there is no mistake
  - Amplicon pipeline: running on GRCh38
  - Sample submission to Genomics
  - 10 projects coming on human cell lines
  - Uploading experimental data: could be automated?
  - Is projectID generated at uploading spreadsheet?
  - Growth: 20 to 50% - filtering out some curves
    - range of colours: for just one category and should be same colour
  - Data table should be ranked by score
  - Better presentation of analysis results: Allele faction and consensus
  - User control: to be added in about 6 month time
  - Adding extra data: post pipeline validation
  - Copy numbers: two variant callers used
  - Add allele frequency plot
  - Clone report: genotype, growth, protein, and all sequence alignments (BAM files with link to IGV)
  - knock-in problem: artificial genome, donor sequence to be store, no cell growth and no

- Future
  - Pool of clones

- ACTIONS:
  - week 16th Oct to load data with Amy
  - 7th November - Amy GitHub training
