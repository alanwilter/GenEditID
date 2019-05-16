import pandas
import logging

from math import ceil, floor

import plotly.graph_objs as go
import plotly.offline as py
from plotly import tools


class Plotter:

    def __init__(self, project_geid):
        self.include_js = False
        self.project_geid = project_geid

    def plate_scoring_plots(self, data_file, plate_scoring_file=None):
            df = pandas.read_csv(data_file)
            # group by plate layout
            dfgroup = df.groupby(['plate'])
            # plot
            # custom colorscale for the score color bar. The scale goes from 0 to 1
            # (0% to 100% of values in the plot)
            color_scale = [[0.0, 'rgb(239, 239, 239)'],  # e.g. this color is applied to the first 10% of values
                           [0.3, 'rgb(200, 223, 247)'],
                           [0.5, 'rgb(237, 119, 90)'],
                           [0.6, 'rgb(237, 195, 90)'],
                           [0.7, 'rgb((237, 229, 90)'],
                           [0.8, 'rgb((181, 221, 24)'],
                           [1, 'rgb(39, 132, 1)']]
            # construct the data dictionary
            dataplot = []
            nloop = 0
            for i, grouped_data in dfgroup:
                nloop += 1  # by default, one scale is created per each trace.
                # I use this counter below to show the first trace's scale and hide the rest
                # create hover text to add to the scatter data structure
                hovertext = []
                for index, rows in grouped_data.iterrows():
                    htext = ["{:s}{:02}".format(rows.row, rows.column), rows.guide, rows.zygosity, 'protein=' + str(round(rows.protein_abundance, 2)), 'score=' + str(rows.score)]
                    htext = ', '.join(htext)  # the hover text needs to be a single string, it can't be a list
                    hovertext.append(htext)
                # create scatter plot data structure
                scatter = go.Scatter(
                    x=grouped_data['column'].tolist(),
                    y=grouped_data['row'].tolist(),
                    name=i,  # i is the plate
                    mode='markers',
                    marker=dict(
                        size='20',  # dot size
                        color=grouped_data['score'].tolist(),  # assign a color based on score
                        showscale=True if nloop == 1 else False,
                        colorscale=color_scale,
                        cmin=7000,  # min value of the colorscale
                        cmax=9000,  # max value of the colorscale
                        colorbar={  # side color bar custom size, fraction of plot size (otherwise it takes all plot height)
                            'lenmode': "fraction",
                            'len': 0.4
                        }
                    ),
                    # hover text
                    text=hovertext,
                    hoverinfo='text'
                )
                dataplot.append(scatter)
            # create plot layout in a grid of two columns and n rows
            # calculate number of subplots (number of plates)
            numberofplates = len(set(df['plate']))
            numberofplotrows = ceil(numberofplates / 2)
            plotheight = numberofplotrows*330  # calculation of total plot height, see comment further down
            # create figure
            subplottitles = [i for i, j in dfgroup]
            figure = tools.make_subplots(rows=numberofplotrows, cols=2, subplot_titles=subplottitles, print_grid=False)
            row_counter = 1
            for i, j in zip(range(0, numberofplates), dataplot):
                subplot_row = floor(row_counter)  # the roughest way to change row every two loops? Oh dear.
                row_counter = row_counter + 0.5
                subplot_column = 1 if 1 & i == 0 else 2
                figure.append_trace(j, subplot_row, subplot_column)
            # update layout construction. In the scatter subplots, each plot has an independent axis
            # (xaxis, xaxis2... xaxisN, and the same for yaxis). Each layout for each subplot
            # must be built independently.
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
            figure.layout.update(dict(autosize=False, width=850, height=plotheight, hovermode='closest', showlegend=False))
            # plot
            output_type = "file"
            if not plate_scoring_file:
                output_type = "div"
                plate_scoring_file = "scores_{}.html".format(self.project_geid)
            return py.plot(figure, filename=plate_scoring_file, auto_open=False, show_link=False,
                           include_plotlyjs=self.include_js, output_type=output_type)


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M')
    try:
        plotter = Plotter('GEP00001')
        plotter.include_js = True
        plotter.plate_scoring_plots("scores.csv", "scores.html")
    except Exception as e:
        logging.exception(e)


if __name__ == '__main__':
    main()
