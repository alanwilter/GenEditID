import plotly.graph_objs as go
import plotly.offline as py

import pandas
import sqlalchemy
import logging

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import VariantResult
from dnascissors.model import Well
from dnascissors.model import WellContent


class Plotter:

    def __init__(self, dbsession, project_geid):
        self.include_js = False
        self.dbsession = dbsession
        self.project_geid = project_geid

    def create_classifier(self, well_content):
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
            classifiers.add(self.create_classifier(well.well_content))
        classifiers = list(classifiers)
        classifiers.sort()
        return classifiers

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
        # Set up colours
        #colours = colorlover.scales['3']['div']['RdYlBu']
        #print(colours)
        #colours = colorlover.interp(colours, len(classifiers))
        #print(colours)
        colour_map = dict()
        colour_index = -1
        for c in classifiers:
            #colour_map[c] = colours[++colour_index]
            colour_map[c] = "blue"
        # Ensure the growths are by increasing time.
        for well in wells:
            well.growths.sort(key=lambda g: g.hours)
        # Need to assemble several plot objects, one for each line.
        # Two loops to order the legend correctly.
        plots = []
        for loop_class in classifiers:
            first = True
            for well in wells:
                classifier = self.create_classifier(well.well_content)
                if classifier == loop_class:
                    plots.append(
                        go.Scatter(
                            mode='lines',
                            line=dict(color=colour_map[classifier]),
                            x=[g.hours for g in well.growths],
                            y=[g.confluence_percentage for g in well.growths],
                            name=classifier,
                            legendgroup=classifier,
                            showlegend=first,
                            hoverinfo='none'
                        )
                    )
                    first = False
        layout = go.Layout(
            title="Cell Growth",
            xaxis=dict(title="Time (h)"),
            yaxis=dict(title="Confluence (%)", range=[0, 100])
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
        #classifiers = self._classifiers_for_wells(wells)
        # Put all the abundances into a list per classifier.
        #classifiers = dict([(c, []) for c in classifiers])
        by_classifier = dict()
        for well in wells:
            c = self.create_classifier(well.well_content)
            l = None
            if c in by_classifier:
                l = by_classifier[c]
            else:
                l = []
                by_classifier[c] = l
            for pa in well.abundances:
                l.append(pa)
        classifiers = list(by_classifier.keys())
        classifiers.sort()
        # Set up colours
        #colours = colorlover.scales['3']['div']['RdYlBu']
        #print(colours)
        #colours = colorlover.interp(colours, len(classifiers))
        #print(colours)
        colour_map = dict()
        colour_index = -1
        for c in classifiers:
            #colour_map[c] = colours[++colour_index]
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
            yaxis=dict(title="Relative protein abundance")
        )
        figure = go.Figure(data=plots, layout=layout)
        output_type = "file"
        if not abundance_file:
            output_type = "div"
            abundance_file = "protein_abundance_{}.html".format(self.project_geid)
        return py.plot(figure, filename=abundance_file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)

    def zygosity_plot(self, zygosity_file=None):
        # we need to filter by gDNA, because there is a sample in project1 (GE-P6B4-G)
        # that is gDNA only (it was sent for sequencing only as gDNA, without a 'fixed cells' counterpart)
        query = self.dbsession.query(Well)\
                              .join(Well.sequencing_library_contents)\
                              .join(Well.experiment_layout)\
                              .join(ExperimentLayout.project)\
                              .filter(Project.geid == self.project_geid)\
                              .filter(SequencingLibraryContent.dna_source != 'gDNA')
        wells = query.all()
        # NB broken
        return None
        if len(wells) == 0:
            return None
        if not wells[0].experiment_layout.project.is_variant_data_available:
            return None
        guides = []
        zygosities = []
        for well in wells:
            layout = well.experiment_layout
            mutation_zygosity = 'wt'
            guide_name = 'none'
            if well.sequencing_library_contents[0].mutation_summaries:
                mutation_zygosity = well.sequencing_library_contents[0].mutation_summaries[0].zygosity
            if well.well_content.guides:
                guide_name = well.well_content.guides[0].name
            zygosities.append(mutation_zygosity)
            guides.append(guide_name)
        # convert 'results' to pandas dataframe and group by 'guides'
        df = pandas.DataFrame({'guides': guides, 'zygosities': zygosities})
        dfgroup = df.groupby(['guides'])

        plots = []
        plots = self.get_percent_plots_by_variant(dfgroup, 'zygosities')
        # order x axis values
        categories = ['wt', 'homo', 'smut', 'dmut', 'iffy']
        layout = go.Layout(
            title='Zygosities',
            xaxis={'categoryorder': 'array', 'categoryarray': categories},
            yaxis={'title': '% of submitted samples per guide'}
        )
        # plot
        figure = go.Figure(data=plots, layout=layout)
        output_type = "file"
        if not zygosity_file:
            output_type = "div"
            zygosity_file = "plot_zygosity_{}.html".format(self.project_geid)
        return py.plot(figure, filename=zygosity_file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)

    def variants_plot(self, variant_type='INDEL', variants_file=None):
        query = self.dbsession.query(VariantResult)\
                              .join(VariantResult.sequencing_library_content)\
                              .join(SequencingLibraryContent.well)\
                              .join(Well.well_content)\
                              .join(Well.experiment_layout)\
                              .join(ExperimentLayout.project)\
                              .filter(Project.geid == self.project_geid)\
                              .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                              .filter(VariantResult.allele_fraction > 0.1)
        variant_results = query.all()
        if len(variant_results) == 0:
            return None
        guides = []
        mutation_types = []
        variant_types = []
        variant_caller_types = []
        for variant_result in variant_results:
            well = variant_result.sequencing_library_content.well
            if not well.experiment_layout.project.is_variant_data_available:
                return None
            guide_name = 'none'
            if well.well_content.guides:
                guide_name = well.well_content.guides[0].name
            guides.append(guide_name)
            mutation_types.append(variant_result.consequence)
            variant_types.append(variant_result.variant_type)
            variant_caller_types.append(variant_result.variant_caller)
        # convert 'variant_results' to pandas dataframe and group by 'variants'
        df = pandas.DataFrame({'variant': variant_types, 'caller': variant_caller_types, 'guide': guides, 'mutation': mutation_types})
        df_group_by_variant = df.groupby(['variant'])

        # This will produce a plot with grouped legends. Annoyingly there is no feature at date 20170710 to add titles to the legend groups.
        # It's been suggested to use layout annotations as a workaround: https://github.com/plotly/plotly.js/issues/689
        # I assume that the top legend group is the first variant caller.
        figure = self.get_percent_plots_by_variant(df_group_by_variant, variant_type, 'mutation')

        output_type = "file"
        if not variants_file:
            output_type = "div"
            variants_file = "variants_{}_{}.html".format(variant_type, self.project_geid)
        return py.plot(figure, filename=variants_file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)

        # the ideal situation would be to plot both datasets side by side
        # pyoff.plot(dict(data = data_indels + data_snvs, layout = layout_indels)) #say we use a different layout for each dataset
        # however the guides appear duplicated. As far as I have seen, you can't have
        # the same legent for two different plots in plotly (ironically you can have it in if you use ggplotly(R:ggplot2)... )

    def get_percent_plots_by_variant(self, group_by_variant, variant_type, grouping_variable):
        plots = []
        marker_symbols = ['circle', 'triangle-up', 'cross', 'hash']
        marker_index = 0
        group_by_variant_type = group_by_variant.get_group(variant_type)
        # calculate percentages per guide in a pandas dataframe
        for variant_caller, group_by_variant_caller in group_by_variant_type.groupby(['caller']):
            for guide_name, group_by_guide in group_by_variant_caller.groupby(['guide']):
                # calculate percentages of grouping_variable and create scatter plot
                group_by_grouping_variable = group_by_guide.groupby([grouping_variable]).size()
                group_by_grouping_variable_percent = group_by_grouping_variable * 100 / group_by_grouping_variable.sum()
                plots.append(
                    go.Scatter(
                         x=group_by_grouping_variable_percent.index.tolist(),
                         y=group_by_grouping_variable_percent.tolist(),
                         name=guide_name,
                         mode='markers',
                         marker={'symbol': marker_symbols[marker_index]},
                         text=variant_caller
                    )
                )
            marker_index += 1
        layout = go.Layout(
                      title='Type of mutation ({}s)'.format(variant_type),
                      yaxis={'title': '% of submitted samples per guide'}
        )
        return go.Figure(data=plots, layout=layout)


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
        plotter.zygosity_plot("zygosity.html")
        plotter.variants_plot("INDEL", "variants_indels.html")
        plotter.variants_plot("SNV", "variants_snvs.html")

    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
