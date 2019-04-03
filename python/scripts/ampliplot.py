import os
import sys
import plotly.offline as py
import plotly.graph_objs as go
import pandas as pd

# output folder name for plots and data
folder_name = 'ampliplots'
if not os.path.exists(folder_name):
    os.mkdir(folder_name)

if not os.path.exists('amplicount_coverage.csv')\
    or not os.path.exists('amplicount.csv')\
    or not os.path.exists('amplicount_config.csv'):
    print('Input files amplicount_config.csv, amplicount_coverage.csv, or amplicount.csv not found, please run amplicount tool first.')
    sys.exit(1)


# Not worse using these two functions yet if there is no characterisation of the variants
def characterise_mutations(mutations):
    if len(mutations) > 0:
        if len(mutations) == 1:
            if mutations[0] > 85:
                return 'homo'
            elif mutations[0] < 85 and mutations[0] > 35:
                return 'smut'
        elif len(mutations) == 2:
            if mutations[0] < 85 and mutations[0] > 35 and mutations[1] < 85 and mutations[1] > 35:
                return 'dmut'
        return 'iffy'


def get_score(mutation_category, consequence):
    score = 0
    if consequence:
        if 'frameshift' in consequence:
            score += 50
    if mutation_category == 'dmut':
        score += 40
    elif mutation_category == 'homo':
        score += 7.5
    elif mutation_category == 'smut':
        score += 2.5
    return score

# Amplicon Read Coverage plot
df_coverage = pd.read_csv('amplicount_coverage.csv')

amplicons_ids = [col for col in df_coverage if col.startswith('chr')]
df_coverage['misc'] = df_coverage['total_reads'] - df_coverage[amplicons_ids].sum(axis=1)

print('Creating Amplicon Read Coverage plot ampliplot.html for {}'.format(amplicons_ids))
data = []
for amplicon_id in amplicons_ids:
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

py.plot({'data': data, 'layout': layout}, filename=os.path.join(folder_name, 'ampliplot.html'), auto_open=False)

# Variant Read Coverage
df_variants = pd.read_csv('amplicount.csv')
df_variants['len'] = df_variants['sequence'].str.len()

df_config = pd.read_csv('amplicount_config.csv')

for amplicon_id in amplicons_ids:
    df_barcodes = df_coverage[['barcode']]
    print('Creating Variant Read Coverage plot {}_ampliplot.html for Amplicon {}'.format(amplicon_id, amplicon_id))
    filter_df_variants = df_variants[(df_variants['amplicon_id'] == amplicon_id) &
                                     (df_variants['len'] > 50)]
    # get reference sequence to mark it as wild-type on plot
    ref_sequence_df = df_variants[(df_variants['amplicon_id'] == amplicon_id) &
                                  (df_variants['seq_type'] == 'REF')]
    if len(ref_sequence_df) > 0:
        ref_sequence = ref_sequence_df.iloc[0]['sequence']
    filter_df_variants = filter_df_variants[['barcode', 'sequence', 'variant_frequency']]
    pivot_df_variants = filter_df_variants.pivot(index='barcode', columns='sequence', values='variant_frequency').reset_index()
    pivot_df_variants.fillna(value=0, inplace=True)
    pivot_df_variants['other'] = 100 - pivot_df_variants.iloc[:, 1:].sum(axis=1)
    pivot_df_variants = pivot_df_variants.merge(df_barcodes, left_on='barcode', right_on='barcode', how='outer')
    pivot_df_variants.fillna(value=0, inplace=True)
    pivot_df_variants.sort_values(by=['barcode'], inplace=True)

    data = []
    i = 0
    with open(os.path.join(folder_name, '{}_ampliplot.fasta'.format(amplicon_id)), 'w') as out:
        for seq in pivot_df_variants.columns.tolist():

            if not seq == 'barcode':
                if seq == 'other':
                    name = seq
                elif seq == ref_sequence:
                    name = 'wild-type'
                    coord = df_config[(df_config['id'] == amplicon_id)].iloc[0]['coord']
                    out.write(">{}_{}_{}\n{}\n".format(amplicon_id, name, coord, ref_sequence))
                else:
                    i += 1
                    name = 'variant_{}'.format(i)
                    out.write(">{}_{}\n{}\n".format(amplicon_id, name, seq))
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
    py.plot({'data': data, 'layout': layout}, filename=os.path.join(folder_name, '{}_ampliplot.html'.format(amplicon_id)), auto_open=False)

# Scores for Amplicons
print('Creating Scores for Amplicons plot ampliscore.html for {}'.format(amplicons_ids))
scores = []
for amplicon_id in amplicons_ids:
    results = []
    for barcode in df_barcodes['barcode'].tolist():
        df_variant_scores = df_variants[(df_variants['amplicon_id'] == amplicon_id) &
                                        (df_variants['barcode'] == barcode) &
                                        (df_variants['len'] > 50) &
                                        (df_variants['seq_type'] != 'REF') &
                                        (df_variants['variant_frequency'] > 20)
                                        ]
        var_type = characterise_mutations(df_variant_scores['variant_frequency'].tolist())
        score = get_score(characterise_mutations(df_variant_scores['variant_frequency'].tolist()), None)
        results.append(score)
        #print(amplicon_id, barcode, len(df_variant_scores), var_type, score)
    scores.append(results)

data = [go.Heatmap(
    x=amplicons_ids,
    y=df_barcodes['barcode'].tolist(),
    z=scores,
    transpose=True,
    xgap=8,
    ygap=1,
    colorscale='Viridis'
)]
layout = go.Layout(
    title='Scores for Amplicons',
    xaxis=dict(
        tickmode='linear'
    )
)
py.plot(go.Figure(data=data, layout=layout), filename=os.path.join(folder_name, 'ampliscore.html'.format(amplicon_id)), auto_open=False)
