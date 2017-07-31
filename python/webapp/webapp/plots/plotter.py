import plotly.graph_objs as go
import plotly.offline as py

import sqlalchemy
import logging

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
                              .filter(Well.growths.any())
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

    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
