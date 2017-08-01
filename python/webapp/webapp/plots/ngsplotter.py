import plotly.graph_objs as go
import plotly.offline as py
from plotly import tools

import pandas
import sqlalchemy
import logging

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import Guide
from dnascissors.model import Target
from dnascissors.model import ExperimentLayout
from dnascissors.model import Project
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import VariantResult
from dnascissors.model import Well
from dnascissors.model import WellContent
from collections import OrderedDict


class NGSPlotter:

    def __init__(self, dbsession, project_geid):
        self.include_js = False
        self.dbsession = dbsession
        self.project_geid = project_geid
        self.legend_groups = set()
        self.allele_fraction_threshold = 0.1
        self.variant_types = ["INDEL", "SNV"]
        self.caller_symbols = {"HaplotypeCaller": "circle", "VarDict": "triangle-up"}
        self.marker_size = 8
        self.guide_color_dict = self._match_colors_to_guidenames()

    def combined_ngs_plot(self, ngs_file=None):
        self.legend_groups.clear()
        titles = []
        for vt in self.variant_types:
            titles.append('Type of mutation ({}s)'.format(vt))
        titles.append("Zygosity")
        titles.append("Indel lenghts")
        # See https://plot.ly/python/subplots/
        self.ngsfigure = tools.make_subplots(rows=4, cols=2,
                                     subplot_titles=titles,
                                     specs=[[{}, {}], [{}, {}], [{'colspan': 2}, None], [{}, {}]],
                                     print_grid=False
                                     )
        self.variants_plot(self.variant_types[0], 1, 1, 1)
        self.variants_plot(self.variant_types[1], 1, 2, 2)
        self.zygosity_plot(2, 1, 3)
        self.indellengths_plot(2, 2, 4)
        self.allelesequences_plot(3, 1, 5)
        output_type = "file"
        if not ngs_file:
            output_type = "div"
            ngs_file = "ngs_{}.html".format(self.project_geid)
        self.ngsfigure.layout.update({
                    'autosize': False,
                    'width': 1000,
                    'height': 1500,
                    'margin': go.Margin(
                        l=50,
                        r=50,
                        b=100,
                        t=100,
                        pad=4)
        })
        return py.plot(self.ngsfigure, filename=ngs_file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)

    def zygosity_plot(self, row_index, column_index, anchor):
        # we need to filter by gDNA, because there is a sample in project1 (GE-P6B4-G)
        # that is gDNA only (it was sent for sequencing only as gDNA, without a 'fixed cells' counterpart)
        query = self.dbsession.query(Well)\
                              .join(Well.sequencing_library_contents)\
                              .join(Well.experiment_layout)\
                              .join(ExperimentLayout.project)\
                              .filter(Project.geid == self.project_geid)\
                              .filter(SequencingLibraryContent.dna_source != 'gDNA')
        wells = query.all()
        if len(wells) == 0:
            return None
        guides = []
        zygosities = []
        for well in wells:
            mutation_zygosity = 'wt'
            guide_name = 'none'
            if well.sequencing_library_contents[0].mutation_summaries:
                mutation_zygosity = well.sequencing_library_contents[0].mutation_summaries[0].zygosity
            if well.well_content.guides:
                guide_name = well.well_content.guides[0].name
            zygosities.append(mutation_zygosity)
            guides.append(guide_name)
        df = pandas.DataFrame({'guides': guides, 'zygosities': zygosities})
        self._get_percent_plots_per_guide(df, 'zygosities', None, row_index, column_index, anchor)

        # order x axis values
        categories = ['wt', 'homo', 'smut', 'dmut', 'iffy']

        self.ngsfigure.layout.update({"xaxis{:d}".format(anchor): dict(categoryorder='array', categoryarray=categories)})
        self.ngsfigure.layout.update({"yaxis{:d}".format(anchor): dict(title='% of submitted samples per guide', range=[0, 100])})

    def variants_plot(self, variant_type, row_index, column_index, anchor):
        query = self.dbsession.query(VariantResult)\
                              .join(VariantResult.sequencing_library_content)\
                              .join(SequencingLibraryContent.well)\
                              .join(Well.well_content)\
                              .join(Well.experiment_layout)\
                              .join(ExperimentLayout.project)\
                              .filter(Project.geid == self.project_geid)\
                              .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                              .filter(VariantResult.allele_fraction > self.allele_fraction_threshold)\
                              .filter(VariantResult.variant_type == variant_type)

        variant_results = query.all()
        if len(variant_results) == 0:
            return None
        guides = []
        mutation_types = []
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
            variant_caller_types.append(variant_result.variant_caller)
        df = pandas.DataFrame({'caller': variant_caller_types, 'guides': guides, 'mutation': mutation_types})
        for variant_caller, group_by_variant_caller in df.groupby(['caller']):
            self._get_percent_plots_per_guide(group_by_variant_caller, 'mutation', variant_caller, row_index, column_index, anchor)
            self.ngsfigure.layout.update({"yaxis{:d}".format(anchor): dict(title='% of submitted samples per guide', range=[0, 100])})

    def indellengths_plot(self, row_index, column_index, anchor):
        query = self.dbsession.query(VariantResult)\
                              .join(VariantResult.sequencing_library_content)\
                              .join(SequencingLibraryContent.well)\
                              .join(Well.experiment_layout)\
                              .join(ExperimentLayout.project)\
                              .filter(Project.geid == self.project_geid)\
                              .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                              .filter(VariantResult.allele_fraction > self.allele_fraction_threshold)\
                              .filter(VariantResult.indel_length.isnot(None))
        variant_results = query.all()
        if len(variant_results) == 0:
            return None
        guides = []
        allindellengths = []
        typecallers = []
        for vr in variant_results:
            well = vr.sequencing_library_content.well
            guide_name = 'none'
            if well.well_content.guides:
                guide_name = well.well_content.guides[0].name
            guides.append(guide_name)
            allindellengths.append(vr.indel_length)
            typecallers.append(vr.variant_caller)
        df = pandas.DataFrame({'caller': typecallers, 'guides': guides, 'indellengths': allindellengths})
        for caller, group_by_variant_caller in df.groupby(['caller']):
            self._get_percent_plots_per_guide(group_by_variant_caller, 'indellengths', caller, row_index, column_index, anchor)
        self.ngsfigure.layout.update({"xaxis{:d}".format(anchor): dict(title='Indel length (bp)')})
        self.ngsfigure.layout.update({"yaxis{:d}".format(anchor): dict(title='% of submitted samples per guide', range=[0, 100])})

    def allelesequences_plot(self, row_index, column_index, anchor):
        query = self.dbsession.query(VariantResult)\
                         .join(VariantResult.sequencing_library_content)\
                         .join(SequencingLibraryContent.well)\
                         .join(Well.well_content)\
                         .join(Well.experiment_layout)\
                         .join(ExperimentLayout.project)\
                         .filter(Project.geid == self.project_geid)\
                         .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                         .filter(VariantResult.allele_fraction > self.allele_fraction_threshold)
        allelesequences = query.all()
        if len(allelesequences) == 0:
            return None
        guide = []
        alleles = []
        typecaller = []
        allele_index5 = []  # this will give the position of '/' in the allele, to then sort the dataframe
        allele_index3 = []
        for i in allelesequences:
            well = i.sequencing_library_content.well
            guide_name = 'none'
            if well.well_content.guides:
                guide_name = well.well_content.guides[0].name
            guide.append(guide_name)
            alleles.append(i.alleles)
            typecaller.append(i.variant_caller)
            allele_split = i.alleles.split('/')
            allele_index5.append(len(allele_split[0]))
            allele_index3.append(len(allele_split[1]))
        # calculate percentages per guide in a pandas dataframe
        # convert 'results' to pandas dataframe and group by 'guides'
        df = pandas.DataFrame({'caller': typecaller, 'guides': guide, 'alleles': alleles, 'index5': allele_index5, 'index3': allele_index3})
        df = df.sort_values(by=['index5', 'index3'])  # calculate percentages of indel lengths and create bar plot 'data' dictionary
        for caller, group_by_variant_caller in df.groupby(['caller']):
            self._get_percent_plots_per_guide(group_by_variant_caller, 'alleles', caller, row_index, column_index, anchor, reverse_axis=False)
        # order the y-axis
        # get the ordered alleles in the y-axis
        allele_list_sorted = list(OrderedDict.fromkeys(df.alleles))
        self.ngsfigure.layout.update({
            "xaxis{:d}".format(anchor): dict(categoryorder='array', categoryarray=allele_list_sorted),
            "yaxis{:d}".format(anchor): dict(title='% of alleles in submitted samples per guide', range=[0, 100]),
            "margin": {'l': 800}
            })

    def _match_colors_to_guidenames(self):
        colour_map_base = ['rgb(0,0,0)', 'rgb(230,159,0)', 'rgb(86,180,233)',
                           'rgb(0,158,115)', 'rgb(240,228,66)', 'rgb(0,114,178)',
                           'rgb(213,94,0)', 'rgb(204,121,167)']
        guidenames = [i.name for i in self.dbsession.query(Guide).join(Target).join(Project).filter(Project.geid == self.project_geid).all()]
        try:
            return dict(zip(guidenames, colour_map_base[0:len(guidenames)]))
        except Exception:
            raise Exception('There are more guides than possible mapping colours for plotting')

    def _get_percent_plots_per_guide(self, df, grouping_variable, variant_caller, row_index, column_index, anchor, reverse_axis = False):
        symbol = self.caller_symbols.get(variant_caller) if variant_caller else 'circle'
        for guide_name, group_by_guide in df.groupby(['guides']):
            colourmap = self.guide_color_dict.get(guide_name)
            group_by_grouping_variable = group_by_guide.groupby([grouping_variable]).size()
            grouped_data_byvar_percent = group_by_grouping_variable*100 / group_by_grouping_variable.sum()
            htext = []
            for length_gdbp in range(len(grouped_data_byvar_percent)):
                if variant_caller:
                    hovertext = [variant_caller] #add other elements to this list to display when hovering
                    hovertext = ' '.join(hovertext)
                    htext.append(hovertext)
                else:
                    hovertext = [] #add other elements to this list to display when hovering
                    hovertext = ' '.join(hovertext)
                    htext.append(hovertext)
            legendgroup = guide_name
            trace = go.Scatter(
                    x = grouped_data_byvar_percent.tolist() if reverse_axis else grouped_data_byvar_percent.index.tolist(),
                    y = grouped_data_byvar_percent.index.tolist() if reverse_axis else grouped_data_byvar_percent.tolist(),
                    name=guide_name,
                    mode='markers',
                    marker=dict(
                        symbol=symbol,
                        size=self.marker_size,
                        color=colourmap
                        ),
                    text=htext,
                    legendgroup=legendgroup,
                    showlegend=legendgroup not in self.legend_groups,
                    xaxis="x{:d}".format(anchor),
                    yaxis="y{:d}".format(anchor))
            print(legendgroup, guide_name, colourmap, variant_caller)
            print(trace)
            self.ngsfigure.append_trace(trace, row_index, column_index)
            self.legend_groups.add(legendgroup)

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M')
    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()
    try:
        plotter = NGSPlotter(session, 'GEP00001')
        plotter.include_js = True
        plotter.combined_ngs_plot("ngs.html")

    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
