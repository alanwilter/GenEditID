import pandas
import sqlalchemy
import logging

import plotly.graph_objs
import plotly.offline
import plotly.tools
import plotly.exceptions

from collections import OrderedDict

from dnascissors.config import cfg
from dnascissors.model import Base
from dnascissors.model import Guide
from dnascissors.model import Target
from dnascissors.model import ExperimentLayout
from dnascissors.model import MutationSummary
from dnascissors.model import Project
from dnascissors.model import SequencingLibraryContent
from dnascissors.model import VariantResult
from dnascissors.model import Well
from dnascissors.model import WellContent


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
        #titles_variant = []
        #for vt in self.variant_types:
        #    titles_variant.append('Consequence of mutation ({}s)'.format(vt))
        #titles.extend(["Indel structure-caller1", "Indel structure-caller2",
        #"Type of mutation", "",
        #'Consequence of mutation (INDELS)', 'Consequence of mutation (SNVs)',
        #"Zygosity", "Indel lenghts",
        #"Allele sequences"])
        # See https://plot.ly/python/subplots/
        titles = ["Type of mutation",
                  "",
                  "Indel structure-caller1",
                  "Indel structure-caller2",
                  "Consequence of mutation (INDELS)",
                  "Consequence of mutation (SNVs)",
                  "Zygosity", "Indel lenghts",
                  "Allele sequences"]
        self.ngsfigure = plotly.tools.make_subplots(rows=5, cols=2,
                                                    subplot_titles=titles,
                                                    specs=[[{}, {}], [{}, {}], [{}, {}], [{},{}], [{'colspan': 2}, None]],
                                                    print_grid=False
                                                    )
        self.typeofmutation_plot(1, 1, 1)
        self.indelstructure_plot(2, [1, 2], [3, 4])  # these plots occupy 2 columns, one anchor per variant caller
        self.variants_plot(self.variant_types[0], 3, 1, 5)
        self.variants_plot(self.variant_types[1], 3, 2, 6)
        self.zygosity_plot(4, 1, 7)
        self.indellengths_plot(4, 2, 8)
        self.allelesequences_plot(5, 1, 9)  # note: anchor 8 gives properties of this plot to indellengths_plot instead

        for anchors in range(1, 10):
            self.ngsfigure.layout.update({"yaxis{:d}".format(anchors): dict(titlefont={'size': 12})})

        output_type = "file"
        if not ngs_file:
            output_type = "div"
            ngs_file = "ngs_{}.html".format(self.project_geid)
        self.ngsfigure.layout.update({
                    'autosize': False,
                    'width': 1000,
                    'height': 1500,
                    'margin': plotly.graph_objs.Margin(
                        l=50,
                        r=50,
                        b=350,
                        t=100,
                        pad=4)
        })
        try:
            return plotly.offline.plot(self.ngsfigure, filename=ngs_file, auto_open=False, show_link=False,
                                       include_plotlyjs=self.include_js, output_type=output_type)
        except plotly.exceptions.PlotlyEmptyDataError:
            return None

    def indelstructure_plot(self, row_index, column_index, anchor, getsample_plot = None):
        """Structure of indels plot.
        It shows the indels and SNVs on a genomic coordinate xaxis, with their size and
        allele frequency. Split by variant caller and showing guide positions.

        """
        query = self.dbsession.query(VariantResult)\
                    .join(VariantResult.sequencing_library_content)\
                    .join(SequencingLibraryContent.well)\
                    .join(SequencingLibraryContent.mutation_summaries)\
                    .join(Well.well_content)\
                    .join(WellContent.guides)\
                    .join(Well.experiment_layout)\
                    .join(ExperimentLayout.project)\
                    .filter(Project.geid == self.project_geid)\
                    .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                    .filter(VariantResult.allele_fraction > self.allele_fraction_threshold)\
                    .filter(MutationSummary.has_off_target == False)
        results = query.all()
        if len(results) == 0:
            return None
        # notes: this query gets only wells with guides, not gDNA, allele fraction > 0.1 and no offtargets

        # get the caller and sequencing_sample_names
        results_callers = list(set(i.variant_caller for i in query.all()))
        results_samples = list(set(i.sequencing_library_content.sequencing_sample_name for i in query.all()))
        # initialise the plot dictionaries
        plotdict = {}
        for caller in results_callers:
            plotcaller = {}
            for sam in results_samples:
                plotcaller.update({sam: []})
            plotdict.update({caller: plotcaller})

        # populate the plot dictionaries on the fly from the database
        plotdata_xaxisrange =[ ]
        guidedict = {}  # list of guides with their coordinates
        shapecolor = {'SNV': 'rgba(0, 0, 0, 0.5)', 'insertion': 'rgba(239, 163, 64, 0.5)', 'deletion': 'rgba(112, 161, 239, 0.5)'}
        anchor_dict = {caller: anc for caller, anc in zip(results_callers, anchor)}
        for i in results:
            well = i.sequencing_library_content.well
            variant_caller = i.variant_caller
            seqsample = i.sequencing_library_content.sequencing_sample_name
            plate = well.experiment_layout.geid
            well_position = "{:s}{:02}".format(well.row, well.column)
            guide_name = well.well_content.guides[0].name
            guide_coordinate = well.well_content.guides[0].amplicon_selections[0].guide_location % 10000 #to get just the last 4 digits, eg. 12345676 > 5676
            xaxisrange = [guide_coordinate - 125, guide_coordinate +125]
            mutation_length = 0 if i.indel_length is None else i.indel_length
            if mutation_length == 0:
                mutation = 'SNV'
            elif mutation_length > 0:
                mutation = 'insertion'
            elif mutation_length < 0:
                mutation = 'deletion'
            mutation_coordinate = i.position
            mutation_coordinate_start = mutation_coordinate if mutation_length >=0 else mutation_coordinate - abs(mutation_length)
            mutation_coordinate_start = mutation_coordinate_start % 10000
            mutation_coordinate_end = mutation_coordinate + mutation_length if mutation_length >= 0 else mutation_coordinate
            mutation_coordinate_end = mutation_coordinate_end % 10000
            allele_fraction = i.allele_fraction
            hovertext = [', '.join([variant_caller, 'well: ' + '-'.join([plate, well_position]),
                                    'sample: ' + seqsample, 'length: '+ str(mutation_length), guide_name])] * 2
            guidedict.update({guide_name: guide_coordinate})
            plotdict[variant_caller][seqsample].append({
                'x': [mutation_coordinate_start, mutation_coordinate_end],
                'y': [allele_fraction]*2,
                'xaxis': "x{:d}".format(anchor_dict[variant_caller]),
                'yaxis': "y{:d}".format(anchor_dict[variant_caller]),
                'name': guide_name,
                'mode': 'lines',
                'legendgroup': guide_name,
                'showlegend': guide_name not in self.legend_groups,
                'line': {
                    'color': shapecolor.get(mutation),
                    'width': 5
                },
                'text': hovertext
            })
            self.legend_groups.add(guide_name)
            plotdata_xaxisrange.append(xaxisrange)

        shapes, annots = [], []
        col_index_dict = {caller: ri for caller, ri in zip(results_callers, column_index)}
        for d_caller in plotdict:
            for d_sample in plotdict[d_caller]:
                for mut in plotdict[d_caller][d_sample]:
                    self.ngsfigure.append_trace(mut, row_index, col_index_dict[d_caller])
            for guidename, guidecoord in zip(guidedict, guidedict.values()):
                shape = {
                    'type': 'line',
                    'xref': "x{:d}".format(anchor_dict[d_caller]),
                    'yref': "y{:d}".format(anchor_dict[d_caller]),
                    'x0': guidecoord, 'y0': 0, 'x1': guidecoord, 'y1': 1,
                    'line': {'color': 'rgb(55, 128, 191)','width': 1, 'dash':'dashdot'},
                }
                annotation = {
                    'text': guidename,
                    'x': guidecoord, 'y':1,
                    'xref': "x{:d}".format(anchor_dict[d_caller]),
                     'yref': "y{:d}".format(anchor_dict[d_caller]),
                    'textangle': 270,
                    'font': {'size':9}
                }
                shapes.append(shape)
                annots.append(annotation)

            self.ngsfigure.layout.update({"yaxis{:d}".format(anchor_dict[d_caller]): {'title': 'mutation structure', 'range': [0, 1.1]},
                                          "xaxis{:d}".format(anchor_dict[d_caller]): {'range': [min(min(plotdata_xaxisrange)), max(max(plotdata_xaxisrange))],
                                                                                      'showgrid': False}
                                          })
        self.ngsfigure.layout.update({'shapes': shapes, 'hovermode': 'closest'})
        self.ngsfigure.layout['annotations'].extend(annots)

        # to get individual sample plots
        fig = plotly.tools.make_subplots(rows=1, cols=2, subplot_titles=results_callers, print_grid=False)
        colindex = 1
        sample_xaxisrange = []
        for d_caller in results_callers:
            for mut in plotdict[d_caller]['GE-P4C2-C']:  # pass getsample_plot here
                fig.append_trace(mut, 1, colindex)
                sample_xaxisrange.append(mut.get('x'))
            colindex += 1

        sample_xaxisranges = [i for x in sample_xaxisrange for i in x]

        for i in fig.layout:
            if i.find('xaxis') == 0:
                fig.layout[i].update({'range': [min(sample_xaxisranges)-10, max(sample_xaxisranges)+10], 'showgrid': False})
            elif i.find('yaxis') == 0:
                fig.layout[i].update({'range': [0, 1.1]})

        fig.layout.update({
            'shapes': shapes,
            'hovermode': 'closest'})
        fig.layout['annotations'].extend(annots)
        #plotly.offline.plot(fig)

    def typeofmutation_plot(self, row_index, column_index, anchor):
        """Type of mutation plot.

        It shows the % of alleles that have a mutation of a certain type.
        It shows allele data from all samples that have a mutation (so wt are excluded,
        and double and single mutants contribute with two and one alleles respectively).
        The types of mutation considered are
        *insertion*, *deletion* and *SNV* (single-nucleotide variant).
        """

        query = self.dbsession.query(VariantResult)\
                              .join(VariantResult.sequencing_library_content)\
                              .join(SequencingLibraryContent.well)\
                              .join(Well.experiment_layout)\
                              .join(ExperimentLayout.project)\
                              .filter(Project.geid == self.project_geid)\
                              .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                              .filter(VariantResult.allele_fraction > self.allele_fraction_threshold)

        varianttype = query.all()
        if len(varianttype) == 0:
            return None
        guides = []
        typeofvariant = []
        variant_caller_types = []
        for i in varianttype:
            guide_name = 'none'
            if i.sequencing_library_content.well.well_content:
                if i.sequencing_library_content.well.well_content.guides:
                    guide_name = i.sequencing_library_content.well.well_content.guides[0].name
            if i.indel_length is None:
                typevar = i.variant_type
            elif i.indel_length > 0:
                typevar = 'insertion'
            elif i.indel_length < 0:
                typevar = 'deletion'
            else:
                typevar = None
            variant_caller_types.append(i.variant_caller)
            guides.append(guide_name)
            typeofvariant.append(typevar)

        df = pandas.DataFrame({'caller': variant_caller_types, 'guides': guides, 'typeofvariant': typeofvariant})
        for variant_caller, group_by_variant_caller in df.groupby(['caller']):
            self._get_percent_plots_per_guide(group_by_variant_caller, 'typeofvariant', variant_caller, row_index, column_index, anchor)
            self.ngsfigure.layout.update({"yaxis{:d}".format(anchor): dict(title='% of submitted samples per guide', range=[0, 100])})

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
            if well.well_content:
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
            if well.well_content:
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
            "xaxis{:d}".format(anchor): dict(categoryorder='array',
                                             categoryarray=allele_list_sorted,
                                             tickfont={'size': 9}),
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

    def _get_percent_plots_per_guide(self, df, grouping_variable, variant_caller, row_index, column_index, anchor, reverse_axis=False):
        symbol = self.caller_symbols.get(variant_caller) if variant_caller else 'circle'
        for guide_name, group_by_guide in df.groupby(['guides']):
            colourmap = self.guide_color_dict.get(guide_name)
            group_by_grouping_variable = group_by_guide.groupby([grouping_variable]).size()
            grouped_data_byvar_percent = group_by_grouping_variable*100 / group_by_grouping_variable.sum()
            htext = []
            for length_gdbp in range(len(grouped_data_byvar_percent)):
                if variant_caller:
                    hovertext = [variant_caller]  # add other elements to this list to display when hovering
                    hovertext = ' '.join(hovertext)
                    htext.append(hovertext)
                else:
                    hovertext = []  # add other elements to this list to display when hovering
                    hovertext = ' '.join(hovertext)
                    htext.append(hovertext)
            legendgroup = guide_name
            trace = plotly.graph_objs.Scatter(
                    x=grouped_data_byvar_percent.tolist() if reverse_axis else grouped_data_byvar_percent.index.tolist(),
                    y=grouped_data_byvar_percent.index.tolist() if reverse_axis else grouped_data_byvar_percent.tolist(),
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
