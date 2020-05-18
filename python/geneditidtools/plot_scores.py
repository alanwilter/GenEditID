import pandas
import glob
import os
import plotly.graph_objs as go
import plotly.offline as py
from plotly import subplots


def main():

    template = pandas.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', '..', 'data', 'templates', 'template_96wellplate.csv'))
    plots_folder = 'geneditid_plots'

    df_all_koscores = pandas.DataFrame()

    for file in glob.glob(os.path.join(plots_folder, 'koscores_*_with_plate_location.csv')):
        amplicon = file.split('koscores_')[1].split('_with_plate_location')[0]
        print(amplicon)
        df_koscores = pandas.read_csv(file)
        df_all_koscores = df_all_koscores.append(df_koscores, ignore_index=True)
        # group by plate layout
        df_koscores_groupby = df_koscores.groupby(['plate_id'])
        # construct the list of data plots
        dataplots = []
        nloop = 0
        for i, grouped_data in df_koscores_groupby:
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
                    size=20,  # dot size
                    color=grouped_data['koscore'].tolist(),  # assign a color based on score
                    showscale=True, 
                    colorscale='Blues',
                    cmin=0,  # min value of the colorscale
                    cmax=1,  # max value of the colorscale
                ),
                # hover text
                text=hovertext,
                hoverinfo='text'
            )
            dataplots.append(scatter)
        # create plot layout in a grid of two columns and n rows
        numberofplates = len(set(df_koscores['plate_id']))
        plotheight = numberofplates*400  # calculation of total plot height, see comment further down
        # create figure
        subplottitles = [i for i, j in df_koscores_groupby]
        figure = subplots.make_subplots(rows=numberofplates, cols=1, subplot_titles=subplottitles, print_grid=False)
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
        py.plot(figure, filename=os.path.join(plots_folder, 'heatmap_{}.html'.format(amplicon)), auto_open=False, show_link=False, include_plotlyjs=True)

    # Protein expression heatmap
    protein_file = 'expression_data.csv'
    if os.path.exists(protein_file):
        df_protein_expression = pandas.read_csv(protein_file)
        df_protein_expression['norm_protein_abundance'] = 1 - (df_protein_expression['protein_abundance'] - df_protein_expression['protein_abundance'].min())/(df_protein_expression['protein_abundance'].max() - df_protein_expression['protein_abundance'].min())
        df_protein_expression = df_protein_expression.merge(df_all_koscores, left_on=['plate_id', 'well', 'sample_id'], right_on=['plate_id', 'well', 'sample_id'], how='left')
        df_protein_expression.fillna(value=0, inplace=True)
        df_protein_expression['combined_score'] = df_protein_expression['norm_protein_abundance']*df_protein_expression['koscore']
        df_protein_expression.to_csv(os.path.join(plots_folder, 'expression_data_normalised_and_combined.csv'), index=False)
        df_protein_expression_groupby = df_protein_expression.groupby(['plate_id'])
        # construct the list of data plots
        dataplots = []
        combined_dataplots = []
        nloop = 0
        for i, grouped_data in df_protein_expression_groupby:
            grouped_data = grouped_data.merge(template, left_on='well', right_on='ref_well', how='right')
            grouped_data = grouped_data[['plate_id', 'sample_id', 'norm_protein_abundance', 'ref_well', 'combined_score']].copy()
            grouped_data.fillna(value={'plate_id': i, 'sample_id': 'no-sample', 'norm_protein_abundance': 0, 'combined_score': 0}, inplace=True)
            grouped_data['row'] = grouped_data['ref_well'].str[1:]
            grouped_data['row'] = grouped_data['row'].astype(int)
            grouped_data['column'] = grouped_data['ref_well'].str[0]
            grouped_data.sort_values(by=['column', 'row'], ascending=[True, True], inplace=True)
            nloop += 1  # by default, one scale is created per each trace.
            # create hover text to add to the scatter data structure
            hovertext = []
            for index, rows in grouped_data.iterrows():
                hovertext.append('{}, {}, protein={}'.format(rows.ref_well, rows.sample_id, rows.norm_protein_abundance))
            # create scatter plot data structure for protein expression
            scatter = go.Scatter(
                x=grouped_data['row'].tolist(),
                y=grouped_data['column'].tolist(),
                name=i,  # i is the plate_id
                mode='markers',
                marker=dict(
                    size='20',  # dot size
                    color=grouped_data['norm_protein_abundance'].tolist(),  # assign a color based on score
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

            # create scatter plot data structure for combined scores
            # create hover text to add to the scatter data structure
            hovertext = []
            for index, rows in grouped_data.iterrows():
                hovertext.append('{}, {}, score={}'.format(rows.ref_well, rows.sample_id, rows.combined_score))
            # create plot
            scatter = go.Scatter(
                x=grouped_data['row'].tolist(),
                y=grouped_data['column'].tolist(),
                name=i,  # i is the plate_id
                mode='markers',
                marker=dict(
                    size='20',  # dot size
                    color=grouped_data['combined_score'].tolist(),  # assign a color based on score
                    showscale=True, # if nloop == 1 else False,
                    colorscale='Reds',
                    reversescale=False,
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
            combined_dataplots.append(scatter)

        # create plot layout in a grid of two columns and n rows
        numberofplates = len(set(df_protein_expression['plate_id']))
        plotheight = numberofplates*400  # calculation of total plot height, see comment further down
        # create figure for protein expression
        subplottitles = [i for i, j in df_protein_expression_groupby]
        figure = subplots.make_subplots(rows=numberofplates, cols=1, subplot_titles=subplottitles, print_grid=False)
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
        py.plot(figure, filename=os.path.join(plots_folder, 'heatmap_protein_expression.html'), auto_open=False, show_link=False, include_plotlyjs=True)

        # create figure for combined plot
        combined_figure = subplots.make_subplots(rows=numberofplates, cols=1, subplot_titles=subplottitles, print_grid=False)
        i = 0
        for dataplot in combined_dataplots:
            i += 1
            combined_figure.append_trace(dataplot, i, 1)
        # update layout construction. In the scatter subplots, each plot has an independent axis
        # (xaxis, xaxis2... xaxisN, and the same for yaxis). Each layout for each subplot must be built independently.
        for i in combined_figure.layout:
            if i.find('xaxis') == 0:
                combined_figure.layout[i].update({'categoryorder': 'array', 'categoryarray': list(range(1, 13))})  # this is to show all xaxis values in the plot, in this order
                combined_figure.layout[i].update({'type': 'category'})
            elif i.find('yaxis') == 0:
                combined_figure.layout[i].update({'categoryorder': 'array', 'categoryarray': ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']})  # this is to show all yaxis values in the plot, in this order
                combined_figure.layout[i].update({'type': 'category'})
                combined_figure.layout[i].update({'type': 'category'})
        # with this layout.update below, width and height are set in the plot. Perhaps you can set them directly on the plotting area on the web page
        # hovermode = closest shows the values for the hover point, otherwise by default ('compare' mode) you only see one coordinate
        combined_figure.layout.update(dict(title='Combined scores', autosize=False, width=600, height=plotheight, hovermode='closest', showlegend=False))
        py.plot(combined_figure, filename=os.path.join(plots_folder, 'heatmap_combined_data.html'), auto_open=False, show_link=False, include_plotlyjs=True)


if __name__ == '__main__':
    main()
