import os
import glob
import pandas as pd

varianid_folder = 'editid_variantid'
df_barcodes = pd.read_csv('barcodes.csv')

df_variants = pd.read_csv(os.path.join(varianid_folder, 'variantid.csv'))
df_impacts = pd.read_csv(os.path.join(varianid_folder, 'impacts.csv'))

df_variants = df_variants.merge(df_barcodes, left_on='sample_id', right_on='sequencing_barcode', how='left')
df_impacts = df_impacts.merge(df_barcodes, left_on='sample_id', right_on='sequencing_barcode', how='left')

df_variants = df_variants[['plate_id', 'well', 'sample_id', 'amplicon_id', 'total_reads', 'amplicon_reads', 'amplicon_filtered_reads', 'amplicon_low_quality_reads', 'amplicon_primer_dimer_reads', 'amplicon_low_abundance_reads', 'variant_reads', 'variant_frequency', 'sequence', 'variant_id', 'variant_type', 'variant_consequence', 'variant_score']]
df_variants.to_csv('editid_variantid/variantid_with_plate_location.csv', index=False)

df_impacts = df_impacts[['plate_id', 'well', 'sample_id', 'amplicon_id', 'impact', 'impact_frequency']]
df_impacts.to_csv('editid_variantid/impacts_with_plate_location.csv', index=False)

for file in glob.glob(os.path.join(varianid_folder, 'koscores_*.csv')):
    if 'with_plate_location' not in file:
        df_koscores = pd.read_csv(file)
        df_koscores = df_koscores.merge(df_barcodes, left_on='sample_id', right_on='sequencing_barcode', how='left')
        df_koscores = df_koscores[['plate_id', 'well', 'sample_id', 'HighImpact', 'MediumImpact', 'LowImpact', 'WildType', 'LowFrequency', 'koscore']]
        output, ext = os.path.splitext(os.path.basename(file))
        df_koscores.to_csv(os.path.join(varianid_folder, '{}_with_plate_location{}'.format(output, ext)), index=False)
