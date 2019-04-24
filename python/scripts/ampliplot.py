import os
import sys
import plotly.offline as py
import plotly.graph_objs as go
import pandas as pd
from Bio import pairwise2
from varcode import Variant
from pyensembl import ensembl_grch38   # pyensembl install --release 95 --species homo_sapiens


# Output folder name for plots and data
folder_name = 'ampliplots'
if not os.path.exists(folder_name):
    os.mkdir(folder_name)

# Checking three output files from amplicount exist
if not os.path.exists('amplicount_coverage.csv')\
    or not os.path.exists('amplicount.csv')\
    or not os.path.exists('amplicount_config.csv'):
    print('Input files amplicount_config.csv, amplicount_coverage.csv, or amplicount.csv not found, please run amplicount tool first.')
    sys.exit(1)

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

# Variant Read Coverage, pairwise alignment to get variant consequences
df_variants = pd.read_csv('amplicount.csv')
df_variants['len'] = df_variants['sequence'].str.len()

df_config = pd.read_csv('amplicount_config.csv')

# Create new output table based on read count outputs
df_variants_output = df_variants[['barcode', 'amplicon_id', 'total_reads', 'amplicon_filtered_reads']]
df_variants_output.rename(columns={'total_reads': 'reads', 'amplicon_filtered_reads': 'reads_filtered'}, inplace=True)
df_variants_output.drop_duplicates(inplace=True)
# wild-type info
df_variants_wt = df_variants[df_variants['seq_type'] == 'REF']
df_variants_wt = df_variants_wt[['barcode', 'amplicon_id', 'variant_reads']]
df_variants_wt.rename(columns={'variant_reads': 'reads_wt'}, inplace=True)
df_variants_output = df_variants_output.merge(df_variants_wt, left_on=['barcode', 'amplicon_id'], right_on=['barcode', 'amplicon_id'], how='left')
df_variants_output.fillna(value=0, inplace=True)
df_variants_output['reads_wt'] = df_variants_output['reads_wt'].astype(int)
# primer-dimer
df_variants_dimer = df_variants[((df_variants['seq_type'] == 'var') & (df_variants['len'] <= 50))]
df_variants_dimer = df_variants_dimer[['barcode', 'amplicon_id', 'variant_reads']]
df_variants_dimer.rename(columns={'variant_reads': 'primer_dimer'}, inplace=True)
sum_df_variants_dimer = df_variants_dimer[['barcode', 'amplicon_id', 'primer_dimer']].groupby(['barcode', 'amplicon_id']).sum()
df_variants_output = df_variants_output.merge(sum_df_variants_dimer, left_on=['barcode', 'amplicon_id'], right_on=['barcode', 'amplicon_id'], how='left')
df_variants_output.fillna(value=0, inplace=True)
df_variants_output['primer_dimer'] = df_variants_output['primer_dimer'].astype(int)
# edited info
df_variants_edited = df_variants[((df_variants['seq_type'] == 'var') & (df_variants['len'] > 50))]
df_variants_edited = df_variants_edited[['barcode', 'amplicon_id', 'sequence', 'variant_reads']]
df_variants_edited.rename(columns={'variant_reads': 'reads_edited'}, inplace=True)
sum_df_variants_edited = df_variants_edited[['barcode', 'amplicon_id', 'reads_edited']].groupby(['barcode', 'amplicon_id']).sum()
count_df_variants_edited = df_variants_edited[['barcode', 'amplicon_id', 'sequence']].groupby(['barcode', 'amplicon_id']).count()
df_variants_output = df_variants_output.merge(count_df_variants_edited, left_on=['barcode', 'amplicon_id'], right_on=['barcode', 'amplicon_id'], how='left')
df_variants_output = df_variants_output.merge(sum_df_variants_edited, left_on=['barcode', 'amplicon_id'], right_on=['barcode', 'amplicon_id'], how='left')
df_variants_output.fillna(value=0, inplace=True)
df_variants_output.rename(columns={'sequence': 'variants'}, inplace=True)
df_variants_output['variants'] = df_variants_output['variants'].astype(int)
df_variants_output['reads_edited'] = df_variants_output['reads_edited'].astype(int)
# define new dataframe for storing sequences with frameshift
df_variants_frameshift = pd.DataFrame(columns=['amplicon_id', 'sequence', 'frameshift', 'deletion', 'insertion'])

#
for amplicon_id in amplicons_ids:
    df_barcodes = df_coverage[['barcode']]
    print('Creating Variant Read Coverage plot {}_ampliplot.html for Amplicon {}'.format(amplicon_id, amplicon_id))
    # select amplicon of interest and filter out primer dimer sequences
    df_variants_filtered = df_variants[(df_variants['amplicon_id'] == amplicon_id) &
                                       (df_variants['len'] > 50)]
    # get reference sequence to mark it as wild-type on plot
    df_ref_sequence = df_variants[(df_variants['amplicon_id'] == amplicon_id) &
                                  (df_variants['seq_type'] == 'REF')]
    if len(df_ref_sequence) > 0:
        ref_sequence = df_ref_sequence.iloc[0]['sequence']
    df_variants_filtered = df_variants_filtered[['barcode', 'sequence', 'variant_frequency']]
    pivot_df_variants_filtered = df_variants_filtered.pivot(index='barcode', columns='sequence', values='variant_frequency').reset_index()
    pivot_df_variants_filtered.fillna(value=0, inplace=True)
    pivot_df_variants_filtered['other'] = 100 - pivot_df_variants_filtered.iloc[:, 1:].sum(axis=1)
    pivot_df_variants_filtered = pivot_df_variants_filtered.merge(df_barcodes, left_on='barcode', right_on='barcode', how='outer')
    pivot_df_variants_filtered.fillna(value=0, inplace=True)
    pivot_df_variants_filtered.sort_values(by=['barcode'], inplace=True)

    data = []
    i = 0
    with open(os.path.join(folder_name, '{}_ampliplot.txt'.format(amplicon_id)), 'w') as out:
        for seq in pivot_df_variants_filtered.columns.tolist():
            if not seq == 'barcode':
                if seq == 'other':
                    name = seq
                elif seq == ref_sequence:
                    name = 'WT'
                    coord = df_config[(df_config['id'] == amplicon_id)].iloc[0]['coord']
                    out.write(">{}_{}_{}\n{}\n".format(amplicon_id, name, coord, ref_sequence))
                else:
                    i += 1
                    #alignments = pairwise2.align.globalms(ref_sequence, seq, 5, -4, -2, -0.5)
                    alignments = pairwise2.align.globalxx(ref_sequence, seq)
                    ibase_start = 0
                    ibase_stop = 0
                    top_effects = set()
                    top_effect_types = set()
                    variant_results = []
                    for ibase in range(1, len(alignments[0][0])-1):
                        top_effect_name = ''
                        top_effect_type = ''
                        if alignments[0][0][ibase] == '-' or alignments[0][1][ibase] == '-':
                            if not alignments[0][0][ibase-1] == '-' and not alignments[0][1][ibase-1] == '-':
                                ibase_start = ibase - 1
                        if ibase_start:
                            if not alignments[0][0][ibase+1] == '-' and not alignments[0][1][ibase+1] == '-':
                                ibase_stop = ibase + 1
                        if ibase_start and ibase_stop:
                            chr, start = amplicon_id.split('_')
                            contig = chr[3:]
                            start = int(start)+ibase_start
                            ref = "{}".format(alignments[0][0][ibase_start:ibase_stop].replace('-', ''))
                            alt = "{}".format(alignments[0][1][ibase_start:ibase_stop].replace('-', ''))
                            try:
                                # Variant consequences using https://github.com/openvax/varcode and https://github.com/openvax/pyensembl
                                var = Variant(contig=contig, start=start, ref=ref, alt=alt, ensembl=ensembl_grch38)
                                top_effect = var.effects().top_priority_effect()
                                top_effect_name = type(top_effect).__name__
                                if top_effect_name == 'FrameShift':
                                    if len(ref) > len(alt):
                                        top_effect_type = 'Deletion'
                                    else:
                                        top_effect_type = 'Insertion'
                            except:
                                top_effect_name = 'Unclassified'
                                pass
                            variant_results.append('{}\t{}\t{}\t{}\t{}'.format(contig, start, ref, alt, top_effect_name))
                            top_effects.add(top_effect_name)
                            top_effect_types.add(top_effect_type)
                            ibase_start = 0
                            ibase_stop = 0
                    name = 'var{}_{}'.format(i, ','.join(top_effects))
                    out.write(">{}_{}\n{}\n".format(amplicon_id, name, seq))
                    out.write(pairwise2.format_alignment(*alignments[0]))
                    out.write('CHROM\tPOS\tREF\tALT\tTOP_EFFECT\n')
                    out.write('{}\n'.format('\n'.join(variant_results)))
                    df_variants_frameshift.loc[len(df_variants_frameshift)] = [amplicon_id, seq, 'FrameShift' in top_effects, 'Deletion' in top_effect_types, 'Insertion' in top_effect_types]
                out.write('\n')
                trace = {
                    'x': pivot_df_variants_filtered[seq],
                    'y': pivot_df_variants_filtered['barcode'],
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

df_variants_frameshift = df_variants_edited.merge(df_variants_frameshift, left_on=['amplicon_id', 'sequence'], right_on=['amplicon_id', 'sequence'], how='left')

df_variants_deletion = df_variants_frameshift[(df_variants_frameshift['deletion'] == True)]
df_variants_deletion.rename(columns={'reads_edited': 'reads_del'}, inplace=True)
groupby_df_variants_deletion = df_variants_deletion[['barcode', 'amplicon_id', 'reads_del']].groupby(['barcode', 'amplicon_id']).sum()

df_variants_insertion = df_variants_frameshift[(df_variants_frameshift['insertion'] == True)]
df_variants_insertion.rename(columns={'reads_edited': 'reads_in'}, inplace=True)
groupby_df_variants_insertion = df_variants_insertion[['barcode', 'amplicon_id', 'reads_in']].groupby(['barcode', 'amplicon_id']).sum()

df_variants_frameshift = df_variants_frameshift[(df_variants_frameshift['frameshift'] == True)]
df_variants_frameshift.rename(columns={'reads_edited': 'reads_frameshift'}, inplace=True)
groupby_df_variants_frameshift = df_variants_frameshift[['barcode', 'amplicon_id', 'reads_frameshift']].groupby(['barcode', 'amplicon_id']).sum()

df_variants_output = df_variants_output.merge(groupby_df_variants_frameshift, left_on=['barcode', 'amplicon_id'], right_on=['barcode', 'amplicon_id'], how='left')
df_variants_output = df_variants_output.merge(groupby_df_variants_deletion, left_on=['barcode', 'amplicon_id'], right_on=['barcode', 'amplicon_id'], how='left')
df_variants_output = df_variants_output.merge(groupby_df_variants_insertion, left_on=['barcode', 'amplicon_id'], right_on=['barcode', 'amplicon_id'], how='left')
df_variants_output.fillna(value=0, inplace=True)
df_variants_output['reads_frameshift'] = df_variants_output['reads_frameshift'].astype(int)
df_variants_output['reads_del'] = df_variants_output['reads_del'].astype(int)
df_variants_output['reads_in'] = df_variants_output['reads_in'].astype(int)
df_variants_output.to_csv(os.path.join(folder_name, 'amplicount_summary.csv'.format(amplicon_id)), index=False)

# Scores for Amplicons
# print('Creating Scores for Amplicons plot ampliscore.html for {}'.format(amplicons_ids))
# scores = []
# for amplicon_id in amplicons_ids:
#     results = []
#     for barcode in df_barcodes['barcode'].tolist():
#         df_variant_scores = df_variants[(df_variants['amplicon_id'] == amplicon_id) &
#                                         (df_variants['barcode'] == barcode) &
#                                         (df_variants['len'] > 50) &
#                                         (df_variants['seq_type'] != 'REF') &
#                                         (df_variants['variant_frequency'] > 20)
#                                         ]
#         var_type = characterise_mutations(df_variant_scores['variant_frequency'].tolist())
#         score = get_score(characterise_mutations(df_variant_scores['variant_frequency'].tolist()), None)
#         results.append(score)
#         #print(amplicon_id, barcode, len(df_variant_scores), var_type, score)
#     scores.append(results)
#
# data = [go.Heatmap(
#     x=amplicons_ids,
#     y=df_barcodes['barcode'].tolist(),
#     z=scores,
#     transpose=True,
#     xgap=8,
#     ygap=1,
#     colorscale='Viridis'
# )]
# layout = go.Layout(
#     title='Scores for Amplicons',
#     xaxis=dict(
#         tickmode='linear'
#     )
# )
# py.plot(go.Figure(data=data, layout=layout), filename=os.path.join(folder_name, 'ampliscore.html'.format(amplicon_id)), auto_open=False)
