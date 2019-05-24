import pandas as pd

df_variants = pd.read_csv('variantid.csv')
df_barcodes = pd.read_csv('barcodes.csv')

df_variants = df_variants.merge(df_barcodes, left_on='sample_id', right_on='sequencing_barcode', how='left')

df_variants = df_variants[['plate_id','well','sample_id','amplicon_id','total_reads','amplicon_reads','amplicon_filtered_reads','amplicon_low_quality_reads','amplicon_primer_dimer_reads','amplicon_low_abundance_reads','variant_reads','variant_frequency','sequence','variant_id','variant_type','variant_consequence','variant_score']]
df_variants.to_csv('variantid_with_plate_location.csv', index=False)
