import plotly.offline as py
import pandas as pd

df_coverage = pd.read_csv('amplicount_coverage.csv')

amlicons_ids = [col for col in df_coverage if col.startswith('chr')]
print(amlicons_ids)
df_coverage['misc'] = df_coverage['total_reads'] - df_coverage[amlicons_ids].sum(axis=1)
print(df_coverage.head())

data = []
for amplicon_id in amlicons_ids:
    trace = {
        'x': df_coverage[amplicon_id],
        'y': df_coverage['barcode'],
        'name': amplicon_id,
        'type': 'bar',
        'orientation': 'h'
    }
    data.append(trace)

misc_trace = {
    'x': df_coverage['misc'],
    'y': df_coverage['barcode'],
    'name': 'misc',
    'type': 'bar',
    'orientation': 'h',
    'marker': {'color': 'rgb(204,204,204)'}
}
data.append(misc_trace)

layout = {'barmode': 'stack',
          'title': 'Amplicon Read Coverage',
          'xaxis': {'title': 'number of reads'},
          'yaxis': {'title': 'samples'}}

py.plot({'data': data, 'layout': layout}, filename='ampliplot')

df_variants = pd.read_csv('amplicount.csv')
df_variants['len'] = df_variants['sequence'].str.len()

print(df_variants.head())

for amplicon_id in amlicons_ids:
    df_barcodes = df_coverage[['barcode']]
    print(amplicon_id)
    filter_df_variants = df_variants[(df_variants['amplicon_id'] == amplicon_id) &\
                                     (df_variants['len'] > 50)]
    # get reference sequence to mark it as wild-type on plot
    ref_sequence_df = df_variants[(df_variants['amplicon_id'] == amplicon_id) &\
                                  (df_variants['seq_type'] == 'REF')]
    if len(ref_sequence_df) > 0:
        ref_sequence = ref_sequence_df.iloc[0]['sequence']
    filter_df_variants = filter_df_variants[['barcode', 'sequence', 'variant_frequency']]
    print(filter_df_variants)
    pivot_df_variants = filter_df_variants.pivot(index='barcode', columns='sequence', values='variant_frequency').reset_index()
    pivot_df_variants.fillna(value=0, inplace=True)
    pivot_df_variants['other'] = 100 - pivot_df_variants.iloc[:, 1:].sum(axis=1)
    pivot_df_variants = pivot_df_variants.merge(df_barcodes, left_on='barcode', right_on='barcode', how='outer')
    pivot_df_variants.fillna(value=0, inplace=True)
    pivot_df_variants.sort_values(by=['barcode'], inplace=True)
    print(pivot_df_variants.head())

    data = []
    i = 0
    for seq in pivot_df_variants.columns.tolist():
        if not seq == 'barcode':
            print(seq)
            if seq == 'other':
                name = seq
            elif seq == ref_sequence:
                name = 'wild-type'
            else:
                i += 1
                name = 'variant_{}'.format(i)
            trace = {
                'x': pivot_df_variants[seq],
                'y': pivot_df_variants['barcode'],
                'name': name,
                'type': 'bar',
                'orientation': 'h'
            }
            data.append(trace)
    layout = {'barmode': 'stack',
              'title': 'Variant Read Coverage for Amplicon {}'.format(amplicon_id),
              'xaxis': {'title': 'frequency of reads'},
              'yaxis': {'title': 'samples'}}
    py.plot({'data': data, 'layout': layout}, filename='{}_ampliplot'.format(amplicon_id))
