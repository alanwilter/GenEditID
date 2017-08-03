import pandas
import sqlalchemy
import logging

from math import ceil, floor

import plotly.graph_objs as go
import plotly.offline as py
from plotly import tools

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import Well
from dnascissors.model import WellContent


class Plotter:

    def __init__(self, dbsession, project_geid):
        self.include_js = False
        self.dbsession = dbsession
        self.project_geid = project_geid

    def growth_plot(self, growth_file=None):
        # See https://stackoverflow.com/questions/21114830/query-to-check-if-size-of-collection-is-0-or-empty-in-sqlalchemy
        # Get all wells for the project where there is growth information.
        query = self.dbsession.query(Well)\
                              .join(Well.well_content)\
                              .join(Well.experiment_layout)\
                              .join(ExperimentLayout.project)\
                              .filter(Project.geid == self.project_geid)\
                              .filter(WellContent.content_type.in_(['sample', 'knock-out', 'wild-type', 'empty-vector', 'normalisation']))\
                              .filter(Well.growths.any())
        wells = query.all()
        if not wells:
            return None
        classifiers = self._classifiers_for_wells(wells)
        # base colour map with 8 colours (OK for colour blind people)
        colour_map = ['rgb(0,0,0)', 'rgb(230,159.0)', 'rgb(86,180,233)',
                      'rgb(0.158,115)', 'rgb(240,228,66)', 'rgb(0,114,178)',
                      'rgb(213,94,0)', 'rgb(204, 121,167)']
        for well in wells:
            well.growths.sort(key=lambda g: g.hours)
        # Need to assemble several plot objects, one for each line.
        # Two loops to order the legend correctly.
        plots = []
        colour_map_index = 0
        for loop_class in classifiers:
            first = True
            for well in wells:
                classifier = self._create_classifier(well.well_content)
                if classifier == loop_class:
                    plots.append(
                        go.Scatter(
                            mode='lines',
                            line=dict(color=colour_map[colour_map_index]),
                            x=[g.hours for g in well.growths],
                            y=[g.confluence_percentage for g in well.growths],
                            name=classifier,
                            legendgroup=classifier,
                            showlegend=first,
                            hoverinfo='none'
                        )
                    )
                    first = False
            colour_map_index += 1
        layout = go.Layout(
            title="Cell Growth",
            xaxis=dict(title="Time (h)"),
            yaxis=dict(title="Confluence (%)", range=[0, 100]),
            autosize=False,
            width=800,
            height=500,
            margin=go.Margin(
                l=50,
                r=50,
                b=100,
                t=100,
                pad=4
            ),
        )
        figure = go.Figure(data=plots, layout=layout)
        output_type = "file"
        if not growth_file:
            output_type = "div"
            growth_file = "cell_growth_{}.html".format(self.project_geid)
        return py.plot(figure, filename=growth_file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)

    def abundance_plot(self, abundance_file=None):
        # Get all wells for the project where there is protein abundance information.
        query = self.dbsession.query(Well)\
                              .join(Well.well_content)\
                              .join(Well.experiment_layout)\
                              .join(ExperimentLayout.project)\
                              .filter(Project.geid == self.project_geid)\
                              .filter(WellContent.content_type.in_(['sample', 'knock-out', 'wild-type', 'normalisation']))\
                              .filter(Well.abundances.any())
        wells = query.all()
        if len(wells) == 0:
            return None
        by_classifier = dict()
        for well in wells:
            c = self._create_classifier(well.well_content)
            l = by_classifier.get(c)
            if not l:
                l = []
                by_classifier[c] = l
            for pa in well.abundances:
                l.append(pa)
        classifiers = list(by_classifier.keys())
        classifiers.sort()
        colour_map = dict()
        for c in classifiers:
            colour_map[c] = "blue"
        # Two loops to order the legend correctly.
        plots = []
        for classifier in classifiers:
            abundances = by_classifier[classifier]
            plots.append(
                go.Scatter(
                    mode='markers',
                    line=dict(color=colour_map[classifier]),
                    x=[classifier] * len(abundances),
                    y=[pa.ratio_800_700 for pa in abundances],
                    name=classifier,
                    legendgroup=classifier,
                    hoverinfo='none'
                )
            )
        layout = go.Layout(
            title="Protein Abundance",
            xaxis=dict(title="Cell Line", ticks=False, fixedrange=True),
            yaxis=dict(title="Relative protein abundance"),
            autosize=False,
            width=800,
            height=500,
            margin=go.Margin(
                l=50,
                r=50,
                b=100,
                t=100,
                pad=4
            ),
        )
        figure = go.Figure(data=plots, layout=layout)
        output_type = "file"
        if not abundance_file:
            output_type = "div"
            abundance_file = "protein_abundance_{}.html".format(self.project_geid)
        return py.plot(figure, filename=abundance_file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)
                       
    def combined_NGSandprotein_plot(self, combinedNGSprot_file=None):
        color_scale = [[0.0, 'rgb(185, 190, 193)'],  # e.g. this color is applied to the first 10% of values
               [0.3, 'rgb(200, 223, 247)'],
               [0.5, 'rgb(237, 119, 90)'],
               [0.6, 'rgb(237, 195, 90)'],
               [0.7, 'rgb((237, 229, 90)'],
               [0.8, 'rgb((181, 221, 24)'],
               [1, 'rgb(39, 132, 1)']]
               
        query = self.dbsession.query(Well)\
                      .join(Well.experiment_layout)\
                      .join(ExperimentLayout.project)\
                      .join(Well.well_content)\
                      .join(WellContent.guides)\
                      .filter(Project.geid == self.project_geid)
        #note I'm getting only the results with guides, ignoring the controls, etc.            
        results = query.all()
        if not results:
            return None
        
        # get all guides in the project
        guidenames = []
        for well in results:
            guidenames.append(well.well_content.guides[0].name)
        
        guidenames = set(guidenames)
        
        # build a list of plot traces
        data = []
        for well in results:
            zygosity = 'notsequenced'
            consequence = 'notsequenced'
            score = 0
            well_abundance_ratio = 0
            
            if well.abundances:
                well_abundance_ratio = well.abundances[0].ratio_800_700
            if well.well_content:
                zygosity = 'wt'
                guide_name = well.well_content.guides[0].name
            if well.sequencing_library_contents:
                # TODO need a loop here to make sure dna_source is fixed cells
                # cannot trust order of items in list returned by sqlalchemy
                # sequencing_library_contents[0].dna_source corresponds to fixed cells,
                # and sequencing_library_contents[1].dna_source to gDNA
                consequence = 'wt'
                if well.sequencing_library_contents[0].mutation_summaries:
                    score = well.sequencing_library_contents[0].mutation_summaries[0].score
                    zygosity = well.sequencing_library_contents[0].mutation_summaries[0].zygosity
                    consequence = well.sequencing_library_contents[0].mutation_summaries[0].consequence
            data_loop = ([
                "{:s}{:02}".format(well.row, well.column),
                well.experiment_layout.geid,
                consequence,
                str(round(well_abundance_ratio, 3)),
                guide_name,
                score,
                zygosity
            ])
            htext = [', '.join([str(i) for i in data_loop])]
            data.append(data_loop + htext)
        
        df = pandas.DataFrame(data)
        column_names = ['wellposition', 'plate', 'consequence', 'protein', 'guide', 'score', 'zygosity', 'hovertext']
        df = df.rename(columns=dict(enumerate(column_names)))
        # group by plate layout
        dfgroup = df.groupby(['guide'])
        
        dataplot = []
        for i, grouped_data in dfgroup:
            scatter = {
                'mode': 'markers',
                'x': grouped_data.consequence.tolist(),
                'y': grouped_data.protein.tolist(),
                'name': i,
                'legendgroup': i,
                'marker': {
                    'color':grouped_data.score.tolist(),  # assign a color based on score
                    'colorscale':color_scale,
                    'cmin':7000,  # min value of the colorscale
                    'cmax':9000,  # max value of the colorscale
                    'colorbar':{  # side color bar custom size, fraction of plot size (otherwise it takes all plot height)
                        'lenmode': "fraction",
                        'len': 0.4
                    }
                 },
                'text': grouped_data.hovertext.tolist()
            }
            dataplot.append(scatter)
        
        categories = sorted(set(df.consequence))
        layout = go.Layout(
            title='Combined NGS and protein',
            xaxis={'categoryorder': 'array', 'categoryarray': categories},
            yaxis={'title': 'Relative protein abundance'},
            hovermode = 'closest'
        )
        
        figure = go.Figure(data = dataplot, layout = layout)
        output_type = "file"
        if not combinedNGSprot_file:
            output_type = "div"
            combinedNGSprot_file = "combinedNGSprot_{}.html".format(self.project_geid)
        return py.plot(figure, filename=combinedNGSprot_file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)
        

    def plate_scoring_plots(self, plate_scoring_file=None):
        results = self.dbsession.query(Well)\
                      .join(Well.experiment_layout)\
                      .join(ExperimentLayout.project)\
                      .filter(Project.geid == self.project_geid)\
                      .all()
        # data collection
        data = []
        for i, well in enumerate(results, start=1):
            zygosity = 'empty-well'
            guide_name = 'no-guide'
            well_content_type = 'empty-well'
            score = 0
            well_abundance_ratio = 0
            if well.abundances:
                well_abundance_ratio = well.abundances[0].ratio_800_700
            if well.well_content:
                well_content_type = 'sample-well'
                zygosity = 'wt'
                if well.well_content.guides:
                    guide_name = well.well_content.guides[0].name
            if well.sequencing_library_contents:
                # TODO need a loop here to make sure dna_source is fixed cells
                # cannot trust order of items in list returned by sqlalchemy
                # sequencing_library_contents[0].dna_source corresponds to fixed cells,
                # and sequencing_library_contents[1].dna_source to gDNA
                if well.sequencing_library_contents[0].mutation_summaries:
                    score = well.sequencing_library_contents[0].mutation_summaries[0].score
                    zygosity = well.sequencing_library_contents[0].mutation_summaries[0].zygosity
            data.append([
                        well.experiment_layout.project.geid,
                        well.experiment_layout.geid,
                        well.row,
                        well.column,
                        guide_name,
                        well_abundance_ratio,
                        zygosity,
                        score,
                        well_content_type
                        ])
        # convert to pandas dataframe and add column names
        df = pandas.DataFrame(data)
        column_names = ['project', 'plate', 'row', 'column', 'guide', 'protein_abundance', 'zygosity', 'score', 'well_content_type']
        df = df.rename(columns=dict(enumerate(column_names)))
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
        # with this layout.update below, width and height are set in the plot.Perhaps you can set them directly on the plotting area on the web page
        # hovermode = closest shows the values for the hover point, otherwise by default ('compare' mode) you only see one coordinate
        figure.layout.update(dict(autosize=False, width=850, height=plotheight, hovermode='closest', showlegend=False))
        # plot
        output_type = "file"
        if not plate_scoring_file:
            output_type = "div"
            plate_scoring_file = "plate_scoring_{}.html".format(self.project_geid)
        return py.plot(figure, filename=plate_scoring_file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)

    def _create_classifier(self, well_content):
        parts = []
        if well_content.clone:
            if well_content.clone.cell_line:
                parts.append(well_content.clone.cell_line.name)
            parts.append(well_content.clone.name)
        if well_content.guides:
            parts.append(",".join(g.name for g in well_content.guides))
        if well_content:
            content = well_content.content_type
            if content and content not in ['empty', 'sample']:
                parts.append(content)
        return " ".join(parts)

    def _classifiers_for_wells(self, wells):
        classifiers = set()
        for well in wells:
            classifiers.add(self._create_classifier(well.well_content))
        classifiers = list(classifiers)
        classifiers.sort()
        return classifiers


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M')
    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()
    try:
        plotter = Plotter(session, 'GEP00001')
        plotter.include_js = True
        plotter.growth_plot("growth.html")
        plotter.abundance_plot("abundance.html")
        plotter.plate_scoring_plots("platescoring.html")
        plotter.combined_NGSandprotein_plot("combined_NGSandprotein_plot.html")

    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
