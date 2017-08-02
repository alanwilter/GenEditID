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
- [ ] add documentation to model.py (e.g. we have a relationship one to many in well to abundances. This is because you can take repeated measurements from the same well)
- [ ] revise relationships model.py (e.g. we have a relationship one to many in sequencing_library_content to mutation_summaries. In this case we can't obtain a second mutation_summaries for the same library, so it's more meaningful if relationship is 1:1).

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
- [x] In the loader, *ProteinAbundanceLoader* class, change file type from '.csv' to '.txt' (or it looks like extension is not used in the loader, only sep?). It's a tab-delimited file anyway, so this way we can keep it consistent with the growth data .txt extension. I am doing the documentation of protein and growth file formats and noting this down. Answer: The extension does not matter for the loader, it only needs to be the same format. We use the csv library to load the file and specify the delimiter to be tab. If we were renaming these files to .txt, we will need to update the loading scripts used to populate the database in `shell/` directory.
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
- RUBEN (1) [x] 96-well plate scatter plot (Ruben - need the scores and the slopes)
- [x] Indelranges (Chandu)
- CHANDU (2) [ ] Type of mutation bar plot (%of samples submitted to NGS vs x =[wt, ins, del, SNV])
- CHANDU (2) [ ] Heatmap [protein, NGS[[frameshift-frame1, frameshift-frame2, noframeshift], [offtarget], [zygosity]], growth slope]
- RUBEN (2) [ ] Growth slopes (we need to calculate the slopes!) - done except with colors
- [ ] Combined plot NGS + protein + slopes, color-coded for frameshifts

- [ ] Not urgent. Current analysis is per single project > update to process multiple projects
- [ ] Not urgent. Current score calculations are based on allele number == 2. However in the PTEN project 3 alleles are knocked, and in future projects multiple positions could be edited simultaneously - need to add an 'allele_number' field to the database and update the code

## pyramid webapp
- **project description**
  - [x] Core Genome Editing comments section: box to add comments about the project, that should be updatable and stored on the database along the project
  - [x] project global overview: tab 'Project' from the layout excel file
  - [x] project detailed overview (on demand, clicking somewhere): tabs 'Project', 'Target', 'Guide' and 'Guidemismatches' from the layout excel file
- **global data exploration**, plots for:
  - RICH [x] Protein abundance
  - RICH [x] Growth curves
  - [ ] Growth slopes
  - NGS exploratory plots
    - ANNE (1) [x] % Zygosities (‘homozygous’, ‘heterozygous’…)
    - [x] % Alleles (e.g. C/CATG, CTAA/C)
    - [ ] % Distances from cut site (to import)
    - [x] % Indel lengths
    - [ ] % Type of mutation (‘wt’, ‘insertion’, ‘deletion’, ‘SNV’)
    - [x] % Type of variant (‘wt’, ‘frameshift’, ‘inframe’...)
    - [ ] INDEL ranges (visual location of alleles)
  - [ ] Combined plots (NGS + protein + growth slopes)
  - [ ] Heatmap
- **sample details**
  - [x] plot: 96-well plates for score-based clone selection
  - RICH [x] Data table per project
  - [ ] Per-sample plot of INDEL ranges
  - LATER [ ] Visualisation of reads (.bam files) in IGV browser (external) (ideally a link to the .bam that opens an IGV server, or an IGV installed locally on the user's computer)
- **user report** for selected samples
  - need a way to select samples and get data table (and plots) only related to these
  - [x] project global overview (tab 'Project' from the layout excel file)
  - [ ] plot INDEL ranges
  - [ ] Combined plot (clone vs controls)
  - Data table
    - [x] per project and plate layout (i.e. GEP00001_01 and no GEP00001_01_incu)
    - [x] output data as soon as we have input (e.g. if there is protein data, show it regardless of presence of growth and NGS data).
    - [ ] Ideally we could click on the 96-well plate and have the clone highlighted in this table, but this is not priority
    - [ ] columns of the data table:
      - [x] Plate
      - [x] Well
      - [x] Sample name
      - [x] Fluidigm barcode
      - [x] score (output NaN if no NGS data available)
      - [ ] Protein abundance (800/700) relative to wild type control (this needs to be calculated. Output NaN if no protein data available)
      - [ ] Protein abundance (800/700) relative to knock-out control (this needs to be calculated)
      - [ ] (Growth slope relative to wild type control (this needs to be calculated. Output NaN if no growth data available))
      - [ ] (Growth slope relative to knock-out control (this needs to be calculated))
      (these below, columns from the variant output. Output NaN if no NGS data available )
      - [x] Variant type/consequence
      - [x] Symbol (Gene ID)
      - [ ] IGV link
      - [x] Allele fraction
      - [x] Alleles
  - [ ] Results Comments (this section is independent to the one described in 'project')
- [ ] generate final user report

## problems to solve
- [x] in plot_typeofvariants.py and plot_distances.py, weird result: for guide STAT3.1 I get a -91 indel_length with haplotypecaller, but according to the excel file that should be STAT3.3 instead (sample GE-P1B5-C).
- [x] is presence of offtargets considered in calculation of zygosities?
- [ ] revise queries, filtering by project might not be working correctly (we can use sample names as tracking system)
- [ ] very large indels detected (90 bp!)
- [ ] issue with *shiny app*: position_dodge messes up hover event data (Ruben is on it, raised issue at github plot.ly)
- CHANDU (1) [x] haplotypecaller is giving two merged values for the same sample onto one row separated by comma
- [ ] unable to fit a model on growth plots to extract curve (Discussion with Dominique on Thu 13 July)
- [ ] make plots and data tables considering that we can have more than one dna_source. As it is currently, e.g. for project GEP00001, well.sequencing_library_contents[0].dna_source is 'fixed cells' and well.sequencing_library_contents[1].dna_source is 'gDNA'. We are selecting only fixed cells for simplification, but we need to show both (the reason to have gDNA and cells was to be able to compare the sequencing results from both). Also, in project one there is a sample (GE-P6B4-G) that is gDNA only (it was sent for sequencing only as gDNA, without a 'fixed cells' counterpart), no fixed cells, so when you do well.sequencing_library_contents[0].dna_source, it results in 'gDNA'!
- [ ] RUBEN. Revise protein data values in plot_96_wells, how the range looks like (we should see edges )

## end of project meeting

Aim: Paper to publish: dry and wet - portable as much as possible: maybe docker


- (1) Bioinfo pipeline:
  - new version of the amplicon sequencing pipeline after fixes on merged data per line
  - branch and modify to our own needs with explanation on how to install from scratch
- (2) Coordinate translation to switch between reference genomes
- (3) Primer design: make it automatic - primer pair if failed nested primers
- Use replicates - update calculation instead of using average
- Question? Several guides in single well: how to you put in the the layout file and load it?
- Document the script methods - 3 to 4 lines of comments in plotter files
- Access control to only genome editing core
  - report with pdf or link to website
  - bam files
- Current webapp - what needs to be finished?
  - [x] Results table need to be updated with the right columns and one line per samples
  - [ ] Add result comments on project and add it on edit project
  - [ ] Finish scoring with protein ratio 800/700 vs KO average
  - [ ] Scale scoring system to 100%
  - [ ] Merge growth and abundance plots into one
  - [x] Check project type: knock-in and knock-out - it is in!
  - [x] Sample tab: update well name on plot to match the one in table
- Update loading data via webapp and ngs results manually
- Load ngs results for knock-in project 3
- Export selection of results not all results
- Add growth slope calculation


- How it went? What we could do better?
  - R/Shiny to Python/Pyramid: is it normal process? nope but lots of learning
  - Agile way of working is great
  - Python, Plotly, database learning
  - Not spend that much time on the shiny
  - 2 times three months deadline has a negative impact on the project flow
  - Allocated working time to work together like two days per week
  - Great to have meeting every morning at 10am, very useful
