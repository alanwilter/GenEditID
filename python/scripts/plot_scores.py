import pandas
import glob
import os
import plotly.graph_objs as go
import plotly.offline as py
from plotly import tools

template = pandas.read_csv('../template_96wellplate.csv')
varianid_folder = 'editid_variantid'

for file in glob.glob(os.path.join(varianid_folder, 'koscores_*_with_plate_location.csv')):
    amplicon = file.split('koscores_')[1].split('_with_plate_location')[0]
    print(amplicon)
    df = pandas.read_csv(file)
    # group by plate layout
    dfgroup = df.groupby(['plate_id'])
    # construct the list of data plots
    dataplots = []
    nloop = 0
    for i, grouped_data in dfgroup:
        grouped_data = grouped_data.merge(template, left_on='well', right_on='ref_well', how='right')
        grouped_data = grouped_data[['plate_id', 'sample_id', 'koscore', 'ref_well']].copy()
        grouped_data.fillna(value={'plate_id': i, 'sample_id': 'no-sample', 'koscore': 0}, inplace=True)
        grouped_data['row'] = grouped_data['ref_well'].str[1:]
        grouped_data['row'] = grouped_data['row'].astype(int)
        grouped_data['column'] = grouped_data['ref_well'].str[0]
        grouped_data.sort_values(by=['column', 'row'], ascending=[True, True], inplace=True)
        nloop += 1  # by default, one scale is created per each trace.
        # create hover text to add to the scatter data structure
        hovertext = []
        for index, rows in grouped_data.iterrows():
            hovertext.append('{}, {}, KOscore={}'.format(rows.ref_well, rows.sample_id, rows.koscore))
        # create scatter plot data structure
        scatter = go.Scatter(
            x=grouped_data['row'].tolist(),
            y=grouped_data['column'].tolist(),
            name=i,  # i is the plate_id
            mode='markers',
            marker=dict(
                size='20',  # dot size
                color=grouped_data['koscore'].tolist(),  # assign a color based on score
                showscale=True, # if nloop == 1 else False,
                colorscale='Blues',
                reversescale=True,
                cmin=0,  # min value of the colorscale
                cmax=100,  # max value of the colorscale
                # colorbar={  # side color bar custom size, fraction of plot size (otherwise it takes all plot height)
                #     'lenmode': "fraction",
                #     'len': 0.6
                # }
            ),
            # hover text
            text=hovertext,
            hoverinfo='text'
        )
        dataplots.append(scatter)
    # create plot layout in a grid of two columns and n rows
    numberofplates = len(set(df['plate_id']))
    plotheight = numberofplates*400  # calculation of total plot height, see comment further down
    # create figure
    subplottitles = [i for i, j in dfgroup]
    figure = tools.make_subplots(rows=numberofplates, cols=1, subplot_titles=subplottitles, print_grid=False)
    i = 0
    for dataplot in dataplots:
        i += 1
        figure.append_trace(dataplot, i, 1)
    # update layout construction. In the scatter subplots, each plot has an independent axis
    # (xaxis, xaxis2... xaxisN, and the same for yaxis). Each layout for each subplot must be built independently.
    for i in figure.layout:
        if i.find('xaxis') == 0:
            figure.layout[i].update({'categoryorder': 'array', 'categoryarray': list(range(1, 13))})  # this is to show all xaxis values in the plot, in this order
            figure.layout[i].update({'type': 'category'})
        elif i.find('yaxis') == 0:
            figure.layout[i].update({'categoryorder': 'array', 'categoryarray': ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']})  # this is to show all yaxis values in the plot, in this order
            figure.layout[i].update({'type': 'category'})
            figure.layout[i].update({'type': 'category'})
    # with this layout.update below, width and height are set in the plot. Perhaps you can set them directly on the plotting area on the web page
    # hovermode = closest shows the values for the hover point, otherwise by default ('compare' mode) you only see one coordinate
    figure.layout.update(dict(title='KO scores for {}'.format(amplicon), autosize=False, width=600, height=plotheight, hovermode='closest', showlegend=False))
    py.plot(figure, filename=os.path.join(varianid_folder, 'heatmap_{}.html'.format(amplicon)), auto_open=False, show_link=False, include_plotlyjs=True)

# Protein expression heatmap
protein_file = 'expression_data_normalised.csv'
if os.path.exists(protein_file):
    df = pandas.read_csv(protein_file)
    dfgroup = df.groupby(['plate_id'])
    # construct the list of data plots
    dataplots = []
    nloop = 0
    for i, grouped_data in dfgroup:
        grouped_data = grouped_data.merge(template, left_on='well', right_on='ref_well', how='right')
        grouped_data = grouped_data[['plate_id', 'sample_id', 'normalised', 'ref_well']].copy()
        grouped_data.fillna(value={'plate_id': i, 'sample_id': 'no-sample', 'normalised': 0}, inplace=True)
        grouped_data['row'] = grouped_data['ref_well'].str[1:]
        grouped_data['row'] = grouped_data['row'].astype(int)
        grouped_data['column'] = grouped_data['ref_well'].str[0]
        grouped_data.sort_values(by=['column', 'row'], ascending=[True, True], inplace=True)
        nloop += 1  # by default, one scale is created per each trace.
        # create hover text to add to the scatter data structure
        hovertext = []
        for index, rows in grouped_data.iterrows():
            hovertext.append('{}, {}, protein={}'.format(rows.ref_well, rows.sample_id, rows.normalised))
        # create scatter plot data structure
        scatter = go.Scatter(
            x=grouped_data['row'].tolist(),
            y=grouped_data['column'].tolist(),
            name=i,  # i is the plate_id
            mode='markers',
            marker=dict(
                size='20',  # dot size
                color=grouped_data['normalised'].tolist(),  # assign a color based on score
                showscale=True, # if nloop == 1 else False,
                colorscale='Greens',
                reversescale=True,
                cmin=0,  # min value of the colorscale
                cmax=1,  # max value of the colorscale
                # colorbar={  # side color bar custom size, fraction of plot size (otherwise it takes all plot height)
                #     'lenmode': "fraction",
                #     'len': 0.6
                # }
            ),
            # hover text
            text=hovertext,
            hoverinfo='text'
        )
        dataplots.append(scatter)
    # create plot layout in a grid of two columns and n rows
    numberofplates = len(set(df['plate_id']))
    plotheight = numberofplates*400  # calculation of total plot height, see comment further down
    # create figure
    subplottitles = [i for i, j in dfgroup]
    figure = tools.make_subplots(rows=numberofplates, cols=1, subplot_titles=subplottitles, print_grid=False)
    i = 0
    for dataplot in dataplots:
        i += 1
        figure.append_trace(dataplot, i, 1)
    # update layout construction. In the scatter subplots, each plot has an independent axis
    # (xaxis, xaxis2... xaxisN, and the same for yaxis). Each layout for each subplot must be built independently.
    for i in figure.layout:
        if i.find('xaxis') == 0:
            figure.layout[i].update({'categoryorder': 'array', 'categoryarray': list(range(1, 13))})  # this is to show all xaxis values in the plot, in this order
            figure.layout[i].update({'type': 'category'})
        elif i.find('yaxis') == 0:
            figure.layout[i].update({'categoryorder': 'array', 'categoryarray': ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']})  # this is to show all yaxis values in the plot, in this order
            figure.layout[i].update({'type': 'category'})
            figure.layout[i].update({'type': 'category'})
    # with this layout.update below, width and height are set in the plot. Perhaps you can set them directly on the plotting area on the web page
    # hovermode = closest shows the values for the hover point, otherwise by default ('compare' mode) you only see one coordinate
    figure.layout.update(dict(title='Protein expression scores', autosize=False, width=600, height=plotheight, hovermode='closest', showlegend=False))
    py.plot(figure, filename=os.path.join(varianid_folder, 'heatmap_protein_expression.html'), auto_open=False, show_link=False, include_plotlyjs=True)


# Combined score heatmap
