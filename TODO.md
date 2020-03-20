# Things done...

- FastQ files added into GDrive to be [downloaded from here](https://drive.google.com/drive/folders/1MN_vzy3hjGOAnycwwtI53nrAuOaB5RJf?usp=sharing) (step 3)
- Access to [Google Drive paper from here](https://drive.google.com/drive/folders/1MQAmhxjuewH2gDoUkzXzz1wmgMK6CV7E?usp=sharing)
- Account setup on [https://framateam.org/geneditid/](https://framateam.org/geneditid/) for team discussion, [please sign-up here!](https://framateam.org/signup_user_complete/?id=q4szea3pitfhpfrmdgaijxsdwo)

# Things in progress...

# Things to do...

- Create proper test data
  - select 3 FastQ files
  - create a simplified submission spreadsheet
  - run all analysis steps on these 3 files only
- Simplify submission spreadsheet and database
- Add plots to the WebApp
  - Add output of the run_ampli_count tool into a new database table for plotting (add class AmpliCountResult in model.py)
  - Update the plotting scripts to call the code from the WebApp to avoid code duplication

- Add validation scripts for the parameters chosen for the alignment `pairwise2.align.globalms(ref_sequence, variant['sequence'], 5, -4, -3, -0.1)`
- Add validation scripts for variant classification `Variant(contig=contig, start=start, ref=ref, alt=alt, ensembl=genome)` `var.effects().top_priority_effect()`
- Modify the analysis steps to be ran from the database instead of flat files
- Drive the analysis from the WebApp
