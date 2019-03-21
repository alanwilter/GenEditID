import plotly.offline as py
import pandas as pd

df = pd.read_csv('amplicount_coverage.csv')

filter_columns = [col for col in df if col.startswith('chr')]
print(filter_columns)
df['misc'] = df['total_reads'] - df[filter_columns].sum(axis=1)
print(df.head())

data = []
for col in filter_columns:
    trace = {
        'x': df[col],
        'y': df['barcode'],
        'name': col,
        'type': 'bar',
        'orientation': 'h'
    }
    data.append(trace)

misc_trace = {
    'x': df['misc'],
    'y': df['barcode'],
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

df = pd.read_csv('amplicount.csv')
print(df.head())
