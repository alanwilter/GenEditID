# Planning meetings

## Data in database
- store NGS results

[Database installation/information](postgres.md)

## What's next? Our TODO list!
# plots / analysis
- [ ] 96-well plate scatter plot (Ruben)
- [ ] Indelranges (Chandu)
- [ ] Type of mutation bar plot (%of samples submitted to NGS vs x =[wt, ins, del, SNV])
- [ ] Heatmap [protein, NGS[[frameshift-frame1, frameshift-frame2, noframeshift], [offtarget], [zygosity]], growth slope]
- [ ] Growth slopes (we need to calculate the slopes!)
- [ ] Combined plot NGS + protein + slopes, color-coded for frameshifts

- [ ] Not urgent. Current analysis is per single project > update to process multiple projects
- [ ] Not urgent. Current score calculations are based on allele number == 2. However in the PTEN project 3 alleles are knocked, and in future projects multiple positions could be edited simultaneously - need to add an 'allele_number' field to the database and update the code
- [ ] Not urgent. currently, one sql query per dataset (NGS, growth, protein). Perhaps one query for all (that allows for missing datasets)

# pipeline
  - [ ] genome coordinates. The user is using hg19 coordinates for primers and guides. We need a script to translate these coordinates to hg18.
  - [ ] pipeline to work with human and mouse genomes
  - [ ] Not urgent. primer design automation and loading data in DB (primer blast)
  
# webapp
- [ ] see docs/webapp-requirements.txt for details on what's essential for the pyramid webapp
- [ ] generate final user report

# excel template 
- [ ] data/templatesYYYYMMDD_GEPXXXXX.xlsx. 
      - we don't need 'barcode_size' in SequencingLibraryContent
      - in tab GuideMismatches: 'regions' would better be 'number_of_mismatches', and 'mismatches' should be 'number_of_offtargets'
      - in Project: change 'institute' to 'affiliation'
      - prepare a knock-in donor table to include in database (not urgent)

# problems to solve
- [ ] in plot_typeofvariants.py and plot_distances.py, weird result: for guide STAT3.1 I get a -91 indel_length with haplotypecaller, but according to the excel file that should be STAT3.3 instead (sample GE-P1B5-C).
- [ ] is presence of offtargets considered in calculation of zygosities?
- [ ] very large indels detected (90 bp!)
- [ ] issue with *shiny app*: position_dodge messes up hover event data (Ruben is on it, raised issue at github plot.ly)
