# Things done...

- FastQ files added into GDrive to be [downloaded from here](https://drive.google.com/drive/folders/1MN_vzy3hjGOAnycwwtI53nrAuOaB5RJf?usp=sharing) (step 3)
- Access to [Google Drive paper from here](https://drive.google.com/drive/folders/1MQAmhxjuewH2gDoUkzXzz1wmgMK6CV7E?usp=sharing)
- Account setup on [https://framateam.org/geneditid/](https://framateam.org/geneditid/) for team discussion, [please sign-up here!](https://framateam.org/signup_user_complete/?id=q4szea3pitfhpfrmdgaijxsdwo)

# Things in progress...

# Things to do...

- Run steps 4 till 6
- Create proper test data
  - simplified submission spreadsheet
  - 3 FastQ files
  - run analysis on these 3 files only for test


- Add plots to the WebApp
- Add output of the run_ampli_count tool into a new database table for plotting in the WebApp (add class AmpliCountResult in model.py)
- Modify the analysis steps to be ran from the database instead of flat files
- Update the plotting scripts to call the code form the WebApp
- Add validation scripts for the parameters chosen for the alignment `pairwise2.align.globalms(ref_sequence, variant['sequence'], 5, -4, -3, -0.1)`
- Add validation scripts for variant classification `Variant(contig=contig, start=start, ref=ref, alt=alt, ensembl=genome)` `var.effects().top_priority_effect()`
- Simplify submission spreadsheet and database
- Drive the analysis from the WebApp
