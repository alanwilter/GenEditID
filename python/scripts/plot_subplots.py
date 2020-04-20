import os
import sys
import math
import pandas
import plotly.offline as py
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def main():
    # Load variants
    variantid_file = os.path.join('data', 'GEP00001', 'editid_variantid', 'variantid.csv')
    df_variants = pandas.read_csv(variantid_file)

    # Amplicon Read Coverage plots
    df_amplicons = df_variants[['sample_id', 'amplicon_id', 'amplicon_reads', 'amplicon_filtered_reads', 'amplicon_low_quality_reads', 'amplicon_primer_dimer_reads', 'amplicon_low_abundance_reads']].copy()
    df_amplicons.drop_duplicates(inplace=True)

    # List of amplicons
    amplicons = df_variants[['amplicon_id']].copy()
    amplicons.drop_duplicates(inplace=True)
    amplicons.reset_index(inplace=True)

    COLORS = {
        'amplicon_filtered_reads': 'rgb(12,100,201)',
        'amplicon_low_quality_reads': 'rgb(204,204,204)',
        'amplicon_primer_dimer_reads': 'rgb(170,170,170)',
        'amplicon_low_abundance_reads': 'rgb(133,133,133)'
    }
    MAX_READS = df_amplicons.loc[df_amplicons['amplicon_reads'].idxmax()]['amplicon_reads']

    titles = []
    for i, amplicon in amplicons.iterrows():
        titles.append('Amplicon Read Coverage for {}'.format(amplicon['amplicon_id']))

    # Initialize figure with subplots
    fig = make_subplots(rows=len(amplicons), cols=1, subplot_titles=titles)

    # Loop over all amplicons
    for i, amplicon in amplicons.iterrows():
        df_coverage = df_amplicons[df_amplicons['amplicon_id'] == amplicon['amplicon_id']].copy()
        df_coverage.sort_values(by=['amplicon_filtered_reads'], inplace=True)

        # Only show legend once
        showlegend=False
        if i == 0:
            showlegend=True

        for name in ['amplicon_filtered_reads', 'amplicon_low_quality_reads', 'amplicon_primer_dimer_reads', 'amplicon_low_abundance_reads']:
            trace = go.Bar(x=df_coverage[name],
                           y=df_coverage['sample_id'],
                           name=' '.join(name.split('_')[1:]),
                           orientation='h',
                           marker={'color': COLORS[name]},
                           showlegend=showlegend,
                          )
            fig.append_trace(trace, i+1, 1)

    layout = {'barmode': 'stack',
              'title': 'Amplicon Read Coverage for {}'.format(amplicon['amplicon_id']),
              'xaxis': {'title': 'number of reads', 'type': 'log', 'range': [0, math.log10(MAX_READS)]},
              'yaxis': {'title': 'samples'}}

    fig.update_layout(barmode='stack', height=800*len(amplicons), width=1200,)
    fig.update_xaxes({'title': 'number of reads', 'type': 'log', 'range': [0, math.log10(MAX_READS)]})
    fig.update_yaxes({'title': 'samples', 'dtick': 1})

    py.plot(fig, filename='coverage.html', auto_open=False)


if __name__ == '__main__':
    main()
