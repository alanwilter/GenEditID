import os
import glob
import math
import pandas
import sqlalchemy
import logging

from Bio import pairwise2
from varcode import Variant
from pyensembl import ensembl_grch38   # pyensembl install --release 95 --species homo_sapiens

import plotly.offline as py
import plotly.graph_objs as go
from plotly import subplots
import plotly.express as px

from geneditid.config import cfg
from geneditid.model import Base
from geneditid.model import LayoutContent
from geneditid.model import Layout
from geneditid.model import Project


class Plotter:

    def __init__(self, dbsession, project_geid):
        self.log = logging.getLogger(__name__)

        self.include_js = False
        self.dbsession = dbsession
        self.project_geid = project_geid

        self.log.info('Plots for project {}'.format(self.project_geid))

        # Folders setup
        self.project_folder = os.path.join(cfg['PROJECTS_FOLDER'], self.project_geid)
        self.plots_folder = os.path.join(self.project_folder, 'geneditid_plots')
        if not os.path.exists(self.plots_folder):
            os.mkdir(self.plots_folder)

        # Consequence configuration
        consequence_config = pandas.read_csv(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'consequences.csv'))

        consequences = consequence_config[['name', 'weight']].copy()
        self.CONSEQUENCE_WEIGHTING = consequences.set_index('name').transpose().to_dict('records')[0]

        categories = consequence_config[['name', 'category']].copy()
        self.CONSEQUENCE_CATEGORIES = categories.set_index('name').transpose().to_dict('records')[0]

        impact = consequence_config[['category', 'weight']].copy()
        impact.drop_duplicates(inplace=True)
        self.IMPACT_WEIGHTING = impact.set_index('category').transpose().to_dict('records')[0]

        # barcodes information
        self.df_barcodes = self.get_barcodes()

        # Load amplicount files (configuration and variants) into pandas dataframe for analysis and ploting
        if self.amplicount_data_exists():
            self.df_config = pandas.read_csv(os.path.join(self.project_folder, 'amplicount_config.csv'))
            self.df_amplicount = pandas.read_csv(os.path.join(self.project_folder, 'amplicount.csv'))
            self.df_variants = self.df_amplicount.copy()
            if self.targeted_search_data_exists():
                self.df_tsearch_config = pandas.read_csv(os.path.join(self.project_folder, 'amplicount_config_tsearch.csv'))
                self.df_tsearch = pandas.read_csv(os.path.join(self.project_folder, 'amplicount_tsearch.csv'))

            if self.df_variants_is_valid():
                # Filter out low-frequency variants
                if len(self.df_variants[(self.df_variants['variant_frequency'] > 5)]) > 0:
                    self.df_variants = self.df_variants[(self.df_variants['variant_frequency'] > 5)]

                # List of sample ids
                self.samples = self.df_variants[['sample_id']].copy()
                self.samples.drop_duplicates(inplace=True)
                self.samples.reset_index(inplace=True)

                # List of amplicon ids
                self.amplicons = self.df_variants[['amplicon_id']].copy()
                self.amplicons.drop_duplicates(inplace=True)
                self.amplicons.reset_index(inplace=True)

                # List of variant ids and sequence
                self.variants = self.df_variants[['amplicon_id', 'sequence']].copy()
                self.variants.drop_duplicates(inplace=True)
                self.variants.reset_index(inplace=True)

                # Classify variants using pairwise alignment
                self.classify_variants()

                # Amplicons dataframe
                self.df_amplicons = self.df_variants[['sample_id', 'amplicon_id', 'amplicon_reads', 'amplicon_filtered_reads', 'amplicon_low_quality_reads', 'amplicon_primer_dimer_reads', 'amplicon_low_abundance_reads']].copy()
                self.df_amplicons.drop_duplicates(inplace=True)

                # Impacts, AllKOScores dataframe
                self.df_impacts = pandas.DataFrame()
                self.df_all_koscores = pandas.DataFrame()


    def get_barcodes(self):
        results = self.dbsession.query(LayoutContent)\
                      .join(LayoutContent.layout)\
                      .join(Layout.project)\
                      .filter(Project.geid == self.project_geid)\
                      .all()
        barcodes = []
        for i, well in enumerate(results, start=1):
            sequencing_barcode = 'no-barcode'
            if well.sequencing_barcode:
                sequencing_barcode = well.sequencing_barcode
            barcodes.append([well.layout.geid,
                             '{}{}'.format(well.row, well.column),
                             sequencing_barcode])

        if barcodes:
            # convert to pandas dataframe and add column names
            df = pandas.DataFrame(barcodes)
            column_names = ['plate_id', 'well', 'sequencing_barcode']
            df = df.rename(columns=dict(enumerate(column_names)))
            df.to_csv(os.path.join(self.project_folder, 'barcodes.csv'), index=False)
            return df


    # Check if amplicount data exists
    def amplicount_data_exists(self):
        if not os.path.exists(os.path.join(self.project_folder, 'amplicount_config.csv')):
            return False
        if not os.path.exists(os.path.join(self.project_folder, 'amplicount.csv')):
            return False
        else:
            if os.stat(os.path.join(self.project_folder, 'amplicount.csv')).st_size == 0:
                return False
        return True


    # Check if desired edits data exists
    def targeted_search_data_exists(self):
        if not os.path.exists(os.path.join(self.project_folder, 'amplicount_config_tsearch.csv')):
            return False
        if not os.path.exists(os.path.join(self.project_folder, 'amplicount_tsearch.csv')):
            return False
        else:
            if os.stat(os.path.join(self.project_folder, 'amplicount_tsearch.csv')).st_size == 0:
                return False
        return True


    # Check Variants dataframe is not empty and has the expected columns
    def df_variants_is_valid(self):
        try:
            if self.df_variants.empty:
                return False
            else:
                for col in ['variant_frequency', 'sample_id', 'sequence', 'amplicon_reads', 'amplicon_filtered_reads', 'amplicon_low_quality_reads', 'amplicon_primer_dimer_reads', 'amplicon_low_abundance_reads']:
                    if col not in self.df_variants.columns:
                        self.log.debug('Missing column {} in amplicount.csv'.format(col))
                        return False
                return True
        except Exception:
            return False


    # Check TargetedSearch dataframe is not empty and has the expected columns
    def df_targeted_search_data_is_valid(self):
        try:
            if self.df_tsearch.empty:
                return False
            else:
                for col in ['sample_id', 'amplicon_id', 'total_reads', 'amplicon_reads', 'amplicon_filtered_reads', 'amplicon_low_quality_reads', 'amplicon_primer_dimer_reads', 'amplicon_low_abundance_reads', 'variant_reads', 'variant_frequency', 'sequence', 'tsearch_id']:
                    if col not in self.df_tsearch.columns:
                        self.log.debug('Missing column {} in amplicount_tsearch.csv'.format(col))
                        return False
                return True
        except Exception:
            return False


    # Get variant classification using https://github.com/openvax/varcode and https://github.com/openvax/pyensembl
    def get_variant_classification(self, contig, start, ref, alt, genome=ensembl_grch38):
        try:
            var = Variant(contig=contig, start=start, ref=ref, alt=alt, ensembl=genome)
            top_effect = var.effects().top_priority_effect()
            consequence = top_effect.__class__.__name__
            weight = self.CONSEQUENCE_WEIGHTING.get(consequence, 0)
        except Exception:
            consequence = 'Unclassified'
            weight = 0
        finally:
            if len(ref) > len(alt):
                return 'Deletion', consequence, weight
            elif len(ref) < len(alt):
                return 'Insertion', consequence, weight
            else:
                return 'Mismatch', consequence, weight


    # Calculate KO score
    def calculate_score(self, row):
        score = 0
        for name in self.IMPACT_WEIGHTING.keys():
            score += row[name]*self.IMPACT_WEIGHTING[name]
        return score/100


    # Pairwise alignment to classify variant
    def classify_variants(self):
        if not os.path.exists(os.path.join(self.plots_folder, 'variantid.csv')):
            with open(os.path.join(self.plots_folder, 'variantid.out'), 'w') as out:
                for i, variant in self.variants.iterrows():
                    # get reference sequence
                    amplicon_id = variant['amplicon_id']
                    df_ref_sequence = self.df_config[(self.df_config['id'] == amplicon_id)]
                    ref_sequence = df_ref_sequence.iloc[0]['amplicon']
                    variant_id = 'var{}'.format(i + 1)
                    data = []
                    top_effect_types = set()
                    top_effect_consequences = set()
                    top_effect_scores = []
                    if variant['sequence'] == ref_sequence:
                        type = 'WildType'
                        top_effect_types.add(type)
                        top_effect_consequences.add(type)
                        top_effect_scores.append(0)
                        coord = self.df_config[(self.df_config['id'] == amplicon_id)].iloc[0]['coord']
                        out.write(">{}_{}_{}_{}\n{}\n".format(amplicon_id, variant_id, type, coord, ref_sequence))
                    else:
                        alignments = pairwise2.align.globalms(ref_sequence, variant['sequence'], 5, -4, -3, -0.1)
                        ibase_start = 0
                        ibase_stop = 0
                        variant_results = []
                        for ibase in range(1, len(alignments[0][0])-1):
                            name, chr, ref_start = amplicon_id.split('_')
                            contig = chr[3:]
                            top_effect_consequence = ''
                            top_effect_type = ''
                            # looking for mismatch
                            if not alignments[0][0][ibase] == alignments[0][1][ibase]:
                                if not alignments[0][0][ibase] == '-' and not alignments[0][1][ibase] == '-':
                                    start = int(ref_start)+ibase
                                    ref = "{}".format(alignments[0][0][ibase])
                                    alt = "{}".format(alignments[0][1][ibase])
                                    top_effect_type, top_effect_consequence, score = self.get_variant_classification(contig, start, ref, alt)
                                    variant_results.append('{}\t{}\t{}\t{}\t{}'.format(contig, start, ref, alt, top_effect_consequence))
                                    top_effect_types.add(top_effect_type)
                                    top_effect_consequences.add(top_effect_consequence)
                                    top_effect_scores.append(score)
                            # looking for insertion and deletion
                            if alignments[0][0][ibase] == '-' or alignments[0][1][ibase] == '-':
                                if not alignments[0][0][ibase-1] == '-' and not alignments[0][1][ibase-1] == '-':
                                    ibase_start = ibase - 1
                            if ibase_start:
                                if not alignments[0][0][ibase+1] == '-' and not alignments[0][1][ibase+1] == '-':
                                    ibase_stop = ibase + 1
                            if ibase_start and ibase_stop:
                                start = int(ref_start)+ibase_start
                                ref = "{}".format(alignments[0][0][ibase_start:ibase_stop].replace('-', ''))
                                alt = "{}".format(alignments[0][1][ibase_start:ibase_stop].replace('-', ''))
                                top_effect_type, top_effect_consequence, score = self.get_variant_classification(contig, start, ref, alt)
                                variant_results.append('{}\t{}\t{}\t{}\t{}'.format(contig, start, ref, alt, top_effect_consequence))
                                top_effect_types.add(top_effect_type)
                                top_effect_consequences.add(top_effect_consequence)
                                top_effect_scores.append(score)
                                ibase_start = 0
                                ibase_stop = 0
                    if len(top_effect_consequences) > 1:
                        if 'FrameShift' in top_effect_consequences:
                            top_effect_consequence = 'ComplexFrameShift'
                            score = self.CONSEQUENCE_WEIGHTING.get(top_effect_consequence, 0)
                        else:
                            top_effect_consequence = 'Complex'
                            score = self.CONSEQUENCE_WEIGHTING.get(top_effect_consequence, 0)
                    else:
                        top_effect_consequence = ''.join(top_effect_consequences)
                        score = top_effect_scores[0]
                    name = '{}_{}'.format(variant_id, ','.join(top_effect_consequences))
                    out.write(">{}_{}\n{}\n".format(amplicon_id, name, variant['sequence']))
                    if variant['sequence'] != ref_sequence:
                        out.write(pairwise2.format_alignment(*alignments[0]))
                        out.write('CHROM\tPOS\tREF\tALT\tTOP_EFFECT\n')
                        out.write('{}\n'.format('\n'.join(variant_results)))
                    self.df_variants.loc[((self.df_variants['amplicon_id'] == amplicon_id) & (self.df_variants['sequence'] == variant['sequence'])), 'variant_id'] = variant_id
                    self.df_variants.loc[((self.df_variants['amplicon_id'] == amplicon_id) & (self.df_variants['sequence'] == variant['sequence'])), 'variant_type'] = '_'.join(top_effect_types)
                    self.df_variants.loc[((self.df_variants['amplicon_id'] == amplicon_id) & (self.df_variants['sequence'] == variant['sequence'])), 'variant_consequence'] = top_effect_consequence
                    self.df_variants.loc[((self.df_variants['amplicon_id'] == amplicon_id) & (self.df_variants['sequence'] == variant['sequence'])), 'variant_score'] = score*self.df_variants['variant_frequency']
                    out.write('\n')

            self.df_variants.drop_duplicates(inplace=True)
            self.df_variants.to_csv(os.path.join(self.plots_folder, 'variantid.csv'), index=False)
        else:
            self.df_variants = pandas.read_csv(os.path.join(self.plots_folder, 'variantid.csv'))


    def config_header(self):
        return list(self.df_config.columns)


    def config_rows(self):
        return self.df_config.to_numpy().tolist()


    def variant_header(self):
        return list(self.df_variants.columns)


    def variant_rows(self):
        return self.df_variants.to_numpy().tolist()


    def amplicount_header(self):
        return list(self.df_amplicount.columns)


    def amplicount_rows(self):
        return self.df_amplicount.to_numpy().tolist()


    def coverage_plot(self, plot_file='coverage.html'):
        if not self.df_variants_is_valid():
            return
        COLORS = {
            'amplicon_filtered_reads': 'rgb(12,100,201)',
            'amplicon_low_quality_reads': 'rgb(204,204,204)',
            'amplicon_primer_dimer_reads': 'rgb(170,170,170)',
            'amplicon_low_abundance_reads': 'rgb(133,133,133)'
        }
        MAX_READS = self.df_amplicons.loc[self.df_amplicons['amplicon_reads'].idxmax()]['amplicon_reads']

        # Plot titles
        titles = []
        for i, amplicon in self.amplicons.iterrows():
            titles.append('Amplicon Read Coverage for {}'.format(amplicon['amplicon_id']))

        # Initialize figure with subplots
        fig = subplots.make_subplots(rows=len(self.amplicons), cols=1, subplot_titles=titles)

        # Loop over all amplicons
        for i, amplicon in self.amplicons.iterrows():
            df_coverage = self.df_amplicons[self.df_amplicons['amplicon_id'] == amplicon['amplicon_id']].copy()
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

        fig.update_layout(barmode='stack', height=1200*len(self.amplicons), width=1200, font=dict(size=8),)
        fig.update_xaxes({'title': 'number of reads', 'type': 'log', 'range': [0, math.log10(MAX_READS)]})
        fig.update_yaxes({'title': 'samples', 'dtick': 1})
        fig.write_html(file=os.path.join(self.plots_folder, plot_file), auto_play=False)
        return fig.to_html(include_plotlyjs=self.include_js, full_html=False)


    def coverage_header(self):
        return list(self.df_amplicons.columns)


    def coverage_rows(self):
        return self.df_amplicons.to_numpy().tolist()


    def variant_impact_plot(self, plot_file='impacts.html'):
        if not self.df_variants_is_valid():
            return

        VCOLORS = {
            'HighImpact': 'rgb(174,19,36)',
            'MediumImpact': 'rgb(206,123,18)',
            'LowImpact': 'rgb(233,185,28)',
            'WildType': 'rgb(250,253,225)',
            'LowFrequency': 'rgb(238,238,238)'
        }

        for consequence in self.CONSEQUENCE_CATEGORIES.keys():
            self.df_variants.loc[(self.df_variants['variant_consequence'] == consequence), 'impact'] = self.CONSEQUENCE_CATEGORIES[consequence]

        self.df_variants.loc[(self.df_variants['variant_consequence'] != 'WildType') & (self.df_variants['variant_score'] == 0), 'impact'] = 'LowImpact'
        self.df_impacts = self.df_variants[['sample_id', 'amplicon_id', 'impact', 'variant_frequency']].copy()
        grouped_impacts = self.df_impacts.groupby(['sample_id', 'amplicon_id', 'impact'])
        self.df_impacts['impact_frequency'] = grouped_impacts.transform('sum')
        self.df_impacts = self.df_impacts.loc[:, ['sample_id', 'amplicon_id', 'impact', 'impact_frequency']]
        self.df_impacts.drop_duplicates(inplace=True)
        self.df_impacts.to_csv(os.path.join(self.plots_folder, 'impacts.csv'), index=False)

        # Plot titles
        titles = []
        for i, amplicon in self.amplicons.iterrows():
            titles.append('Variant Impact Frequency for {}'.format(amplicon['amplicon_id']))

        # Initialize figure with subplots
        fig = subplots.make_subplots(rows=len(self.amplicons), cols=1, subplot_titles=titles)

        # list of koscore dataframe for all amplicons
        all_df_koscores = []

        # Loop over all amplicons
        for i, amplicon in self.amplicons.iterrows():
            df_impacts_per_amplicon = self.df_impacts[self.df_impacts['amplicon_id'] == amplicon['amplicon_id']]
            pivot_df_impacts_per_amplicon = df_impacts_per_amplicon.pivot(index='sample_id', columns='impact', values='impact_frequency').reset_index()
            pivot_df_impacts_per_amplicon.fillna(value=0, inplace=True)
            pivot_df_impacts_per_amplicon['LowFrequency'] = 100 - pivot_df_impacts_per_amplicon.iloc[:, 1:].sum(axis=1)
            for name in self.IMPACT_WEIGHTING.keys():
                if name not in pivot_df_impacts_per_amplicon.columns:
                    pivot_df_impacts_per_amplicon[name] = 0
            pivot_df_impacts_per_amplicon['koscore'] = pivot_df_impacts_per_amplicon.apply(self.calculate_score, axis=1)
            pivot_df_impacts_per_amplicon.sort_values(by=['koscore'], ascending=[True], inplace=True)
            # add barcodes
            df_koscores = pivot_df_impacts_per_amplicon.copy()
            df_koscores = df_koscores.merge(self.df_barcodes, left_on='sample_id', right_on='sequencing_barcode', how='left')
            df_koscores = df_koscores[['plate_id', 'well', 'sample_id', 'HighImpact', 'MediumImpact', 'LowImpact', 'WildType', 'LowFrequency', 'koscore']]
            df_koscores.insert(0, 'amplicon_id', amplicon['amplicon_id'])
            all_df_koscores.append(df_koscores)

            # Only show legend once
            showlegend=False
            if i == 0:
                showlegend=True

            for name in ['HighImpact', 'MediumImpact', 'LowImpact', 'WildType', 'LowFrequency']:
                trace = go.Bar(x=pivot_df_impacts_per_amplicon[name],
                               y=pivot_df_impacts_per_amplicon['sample_id'],
                               name=name,
                               orientation='h',
                               marker={'color': VCOLORS[name]},
                               showlegend=showlegend,
                              )
                fig.append_trace(trace, i+1, 1)

        fig.update_layout(barmode='stack', height=1200*len(self.amplicons), width=1200, font=dict(size=8),)
        fig.update_xaxes({'title': 'frequency'})
        fig.update_yaxes({'title': 'samples', 'dtick': 1})
        fig.write_html(file=os.path.join(self.plots_folder, plot_file), auto_play=False)
        self.df_all_koscores = pandas.concat(all_df_koscores, ignore_index=True)
        self.df_all_koscores.to_csv(os.path.join(self.plots_folder, 'koscores.csv'), index=False)
        return fig.to_html(include_plotlyjs=self.include_js, full_html=False)


    def impact_header(self):
        return list(self.df_impacts.columns)


    def impact_rows(self):
        return self.df_impacts.to_numpy().tolist()


    def heatmap_plot(self, plot_file='koscores.html'):
        if not self.df_variants_is_valid():
            return
        if self.df_all_koscores.empty:
            return

        template = pandas.read_csv(os.path.join('data', 'templates', 'template_96wellplate.csv'))

        # plates
        plates = list(set(self.df_all_koscores['plate_id']))
        plates.sort()

        # amplicons
        amplicons = list(set(self.df_all_koscores['amplicon_id']))
        amplicons.sort()

        # calculate total plot size
        plotheight = len(plates)*400
        plotwidth = len(amplicons)*600

        # Plot titles
        titles = []
        for plate in plates:
            for amplicon in amplicons:
                titles.append('{} on {}'.format(amplicon, plate))

        # Initialize figure with subplots
        fig = subplots.make_subplots(rows=len(plates), cols=len(amplicons), subplot_titles=titles, print_grid=False)

        # Loop over all amplicons and plates
        for i, amplicon in enumerate(amplicons, start=1):
            for j, plate in enumerate(plates, start=1):
                df_koscores_on_plate = self.df_all_koscores[(self.df_all_koscores['amplicon_id'] == amplicon) & (self.df_all_koscores['plate_id'] == plate)]
                df_koscores_on_plate = df_koscores_on_plate.merge(template, left_on='well', right_on='ref_well', how='right')
                df_koscores_on_plate = df_koscores_on_plate[['plate_id', 'sample_id', 'koscore', 'ref_well']].copy()
                df_koscores_on_plate.fillna(value={'plate_id': plate, 'sample_id': 'no-sample', 'koscore': 0}, inplace=True)
                df_koscores_on_plate['row'] = df_koscores_on_plate['ref_well'].str[1:]
                df_koscores_on_plate['row'] = df_koscores_on_plate['row'].astype(int)
                df_koscores_on_plate['column'] = df_koscores_on_plate['ref_well'].str[0]
                df_koscores_on_plate.sort_values(by=['column', 'row'], ascending=[True, True], inplace=True)
                # create hover text to add to the scatter data structure
                hovertext = []
                for index, rows in df_koscores_on_plate.iterrows():
                    hovertext.append('{}, {}, KOscore={}'.format(rows.ref_well, rows.sample_id, rows.koscore))
                # create scatter plot data structure
                trace = go.Scatter(
                    x=df_koscores_on_plate['row'].tolist(),
                    y=df_koscores_on_plate['column'].tolist(),
                    name=plate,
                    mode='markers',
                    marker=dict(
                        size=20,  # dot size
                        color=df_koscores_on_plate['koscore'].tolist(),  # assign a color based on score
                        showscale=True,
                        colorscale='Blues',
                        cmin=0,  # min value of the colorscale
                        cmax=1,  # max value of the colorscale
                    ),
                    # hover text
                    text=hovertext,
                    hoverinfo='text'
                )
                fig.append_trace(trace, j, i)
                # Update plot layout construction to be in a grid of two columns and n rows
                # In the scatter subplots, each plot has an independent axis
                # (xaxis, xaxis2... xaxisN, and the same for yaxis). Each layout for each subplot must be built independently.
                for l in fig.layout:
                    if l.find('xaxis') == 0:
                        fig.layout[l].update({'categoryorder': 'array', 'categoryarray': list(range(1, 13))})  # this is to show all xaxis values in the plot, in this order
                        fig.layout[l].update({'type': 'category'})
                    elif l.find('yaxis') == 0:
                        fig.layout[l].update({'categoryorder': 'array', 'categoryarray': ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']})  # this is to show all yaxis values in the plot, in this order
                        fig.layout[l].update({'type': 'category'})
                        fig.layout[l].update({'type': 'category'})
                # with this layout.update below, width and height are set in the plot
                # hovermode = closest shows the values for the hover point, otherwise by default ('compare' mode) you only see one coordinate
        fig.layout.update(dict(autosize=False, width=plotwidth, height=plotheight, hovermode='closest', showlegend=False))
        fig.write_html(file=os.path.join(self.plots_folder, plot_file), auto_play=False)
        return fig.to_html(include_plotlyjs=self.include_js, full_html=False)


    def koscores_header(self):
        return list(self.df_all_koscores.columns)


    def koscores_rows(self):
        return self.df_all_koscores.to_numpy().tolist()


    def targeted_search_plot(self, plot_file='targeted_search.html'):
        if not self.targeted_search_data_exists():
            return
        if not self.df_targeted_search_data_is_valid():
            return

        MAX_READS = self.df_tsearch.loc[self.df_tsearch['amplicon_filtered_reads'].idxmax()]['amplicon_filtered_reads']

        # List of amplicon ids
        amplicons = self.df_tsearch[['amplicon_id']].copy()
        amplicons.drop_duplicates(inplace=True)
        amplicons.reset_index(inplace=True)

        # List of targeted search ids
        tsearches = self.df_tsearch[['tsearch_id']].copy()
        tsearches.drop_duplicates(inplace=True)
        tsearches.reset_index(inplace=True)

        df_tcoverage_frequency = self.df_tsearch[['sample_id', 'amplicon_id', 'variant_frequency', 'tsearch_id']].copy()
        gby_tcoverage_frequency = df_tcoverage_frequency.groupby(
                                        ['sample_id', 'amplicon_id', 'tsearch_id'],
                                        as_index=False).agg({'variant_frequency' : sum })
        gby_tcoverage_frequency = gby_tcoverage_frequency.loc[:, ['sample_id', 'amplicon_id', 'tsearch_id', 'variant_frequency']]

        # Plot titles
        titles = []
        for i, amplicon in amplicons.iterrows():
            titles.append('Targeted Search for Amplicon {}'.format(amplicon['amplicon_id']))

        # Initialize figure with subplots
        fig = subplots.make_subplots(rows=len(amplicons), cols=1, subplot_titles=titles)

        # Loop over all amplicons
        for i, amplicon in amplicons.iterrows():
            gby_tcoverage_per_amplicon = gby_tcoverage_frequency[gby_tcoverage_frequency['amplicon_id'] == amplicon['amplicon_id']].copy()
            pivot_gby_tcoverage_per_amplicon = gby_tcoverage_per_amplicon.pivot(index='sample_id', columns='tsearch_id', values='variant_frequency').reset_index()
            pivot_gby_tcoverage_per_amplicon.fillna(value=0, inplace=True)
            pivot_gby_tcoverage_per_amplicon['Others'] = 100 - pivot_gby_tcoverage_per_amplicon.iloc[:, 1:].sum(axis=1)
            pivot_gby_tcoverage_per_amplicon.sort_values(by=['Others'], ascending=[False], inplace=True)
            for name in tsearches['tsearch_id'].tolist():
                if name not in pivot_gby_tcoverage_per_amplicon.columns:
                    pivot_gby_tcoverage_per_amplicon[name] = 0

            # Only show legend once
            showlegend=False
            if i == 0:
                showlegend=True
            j = 0
            for name in tsearches['tsearch_id'].tolist():
                trace = go.Bar(x=pivot_gby_tcoverage_per_amplicon[name],
                               y=pivot_gby_tcoverage_per_amplicon['sample_id'],
                               name=name,
                               orientation='h',
                               marker={'color': px.colors.qualitative.Antique[j%10]},
                               showlegend=showlegend,
                              )
                fig.append_trace(trace, i+1, 1)
                j += 1
            trace = go.Bar(x=pivot_gby_tcoverage_per_amplicon['Others'],
                           y=pivot_gby_tcoverage_per_amplicon['sample_id'],
                           name='Others',
                           orientation='h',
                           marker={'color': 'rgb(204,204,204)'},  # grey for Others
                           showlegend=showlegend,
                          )
            fig.append_trace(trace, i+1, 1)

        fig.update_layout(barmode='stack', height=1200*len(amplicons), width=1200, font=dict(size=8),)
        fig.update_xaxes({'title': 'frequency of filtered reads'})
        fig.update_yaxes({'title': 'samples', 'dtick': 1})
        fig.write_html(file=os.path.join(self.plots_folder, plot_file), auto_play=False)
        return fig.to_html(include_plotlyjs=self.include_js, full_html=False)


    def tsearch_config_header(self):
        return list(self.df_tsearch_config.columns)


    def tsearch_config_rows(self):
        return self.df_tsearch_config.to_numpy().tolist()


    def tsearch_header(self):
        return list(self.df_tsearch.columns)


    def tsearch_rows(self):
        return self.df_tsearch.to_numpy().tolist()


        # # Protein expression heatmap
        # protein_file = os.path.join('data', 'GEP00001', 'expression_data.csv')
        # if os.path.exists(protein_file):
        #     df_protein_expression = pandas.read_csv(protein_file)
        #     df_protein_expression['norm_protein_abundance'] = 1 - (df_protein_expression['protein_abundance'] - df_protein_expression['protein_abundance'].min())/(df_protein_expression['protein_abundance'].max() - df_protein_expression['protein_abundance'].min())
        #     df_protein_expression = df_protein_expression.merge(df_all_koscores, left_on=['plate_id', 'well', 'sample_id'], right_on=['plate_id', 'well', 'sample_id'], how='left')
        #     df_protein_expression.fillna(value=0, inplace=True)
        #     df_protein_expression['combined_score'] = df_protein_expression['norm_protein_abundance']*df_protein_expression['koscore']
        #     df_protein_expression.to_csv(os.path.join(self.plots_folder, 'expression_data_normalised_and_combined.csv'), index=False)
        #     df_protein_expression_groupby = df_protein_expression.groupby(['plate_id'])
        #     # construct the list of data plots
        #     dataplots = []
        #     combined_dataplots = []
        #     nloop = 0
        #     for i, grouped_data in df_protein_expression_groupby:
        #         grouped_data = grouped_data.merge(template, left_on='well', right_on='ref_well', how='right')
        #         grouped_data = grouped_data[['plate_id', 'sample_id', 'norm_protein_abundance', 'ref_well', 'combined_score']].copy()
        #         grouped_data.fillna(value={'plate_id': i, 'sample_id': 'no-sample', 'norm_protein_abundance': 0, 'combined_score': 0}, inplace=True)
        #         grouped_data['row'] = grouped_data['ref_well'].str[1:]
        #         grouped_data['row'] = grouped_data['row'].astype(int)
        #         grouped_data['column'] = grouped_data['ref_well'].str[0]
        #         grouped_data.sort_values(by=['column', 'row'], ascending=[True, True], inplace=True)
        #         nloop += 1  # by default, one scale is created per each trace.
        #         # create hover text to add to the scatter data structure
        #         hovertext = []
        #         for index, rows in grouped_data.iterrows():
        #             hovertext.append('{}, {}, protein={}'.format(rows.ref_well, rows.sample_id, rows.norm_protein_abundance))
        #         # create scatter plot data structure for protein expression
        #         scatter = go.Scatter(
        #             x=grouped_data['row'].tolist(),
        #             y=grouped_data['column'].tolist(),
        #             name=i,  # i is the plate_id
        #             mode='markers',
        #             marker=dict(
        #                 size=20,  # dot size
        #                 color=grouped_data['norm_protein_abundance'].tolist(),  # assign a color based on score
        #                 showscale=True, # if nloop == 1 else False,
        #                 colorscale='Greens',
        #                 cmin=0,  # min value of the colorscale
        #                 cmax=1,  # max value of the colorscale
        #             ),
        #             # hover text
        #             text=hovertext,
        #             hoverinfo='text'
        #         )
        #         dataplots.append(scatter)
        #
        #         # create scatter plot data structure for combined scores
        #         # create hover text to add to the scatter data structure
        #         hovertext = []
        #         for index, rows in grouped_data.iterrows():
        #             hovertext.append('{}, {}, score={}'.format(rows.ref_well, rows.sample_id, rows.combined_score))
        #         # create plot
        #         scatter = go.Scatter(
        #             x=grouped_data['row'].tolist(),
        #             y=grouped_data['column'].tolist(),
        #             name=i,  # i is the plate_id
        #             mode='markers',
        #             marker=dict(
        #                 size=20,  # dot size
        #                 color=grouped_data['combined_score'].tolist(),  # assign a color based on score
        #                 showscale=True, # if nloop == 1 else False,
        #                 colorscale='Reds',
        #                 cmin=0,  # min value of the colorscale
        #                 cmax=1,  # max value of the colorscale
        #             ),
        #             # hover text
        #             text=hovertext,
        #             hoverinfo='text'
        #         )
        #         combined_dataplots.append(scatter)
        #
        #     # create plot layout in a grid of two columns and n rows
        #     numberofplates = len(set(df_protein_expression['plate_id']))
        #     plotheight = numberofplates*400  # calculation of total plot height, see comment further down
        #     # create figure for protein expression
        #     subplottitles = [i for i, j in df_protein_expression_groupby]
        #     figure = subplots.make_subplots(rows=numberofplates, cols=1, subplot_titles=subplottitles, print_grid=False)
        #     i = 0
        #     for dataplot in dataplots:
        #         i += 1
        #         figure.append_trace(dataplot, i, 1)
        #     # update layout construction. In the scatter subplots, each plot has an independent axis
        #     # (xaxis, xaxis2... xaxisN, and the same for yaxis). Each layout for each subplot must be built independently.
        #     for i in figure.layout:
        #         if i.find('xaxis') == 0:
        #             figure.layout[i].update({'categoryorder': 'array', 'categoryarray': list(range(1, 13))})  # this is to show all xaxis values in the plot, in this order
        #             figure.layout[i].update({'type': 'category'})
        #         elif i.find('yaxis') == 0:
        #             figure.layout[i].update({'categoryorder': 'array', 'categoryarray': ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']})  # this is to show all yaxis values in the plot, in this order
        #             figure.layout[i].update({'type': 'category'})
        #             figure.layout[i].update({'type': 'category'})
        #     # with this layout.update below, width and height are set in the plot. Perhaps you can set them directly on the plotting area on the web page
        #     # hovermode = closest shows the values for the hover point, otherwise by default ('compare' mode) you only see one coordinate
        #     figure.layout.update(dict(title='Protein expression scores', autosize=False, width=600, height=plotheight, hovermode='closest', showlegend=False))
        #     py.plot(figure, filename=os.path.join(self.plots_folder, 'heatmap_protein_expression.html'), auto_open=False, show_link=False, include_plotlyjs=True)
        #
        #     # create figure for combined plot
        #     combined_figure = subplots.make_subplots(rows=numberofplates, cols=1, subplot_titles=subplottitles, print_grid=False)
        #     i = 0
        #     for dataplot in combined_dataplots:
        #         i += 1
        #         combined_figure.append_trace(dataplot, i, 1)
        #     # update layout construction. In the scatter subplots, each plot has an independent axis
        #     # (xaxis, xaxis2... xaxisN, and the same for yaxis). Each layout for each subplot must be built independently.
        #     for i in combined_figure.layout:
        #         if i.find('xaxis') == 0:
        #             combined_figure.layout[i].update({'categoryorder': 'array', 'categoryarray': list(range(1, 13))})  # this is to show all xaxis values in the plot, in this order
        #             combined_figure.layout[i].update({'type': 'category'})
        #         elif i.find('yaxis') == 0:
        #             combined_figure.layout[i].update({'categoryorder': 'array', 'categoryarray': ['H', 'G', 'F', 'E', 'D', 'C', 'B', 'A']})  # this is to show all yaxis values in the plot, in this order
        #             combined_figure.layout[i].update({'type': 'category'})
        #             combined_figure.layout[i].update({'type': 'category'})
        #     # with this layout.update below, width and height are set in the plot. Perhaps you can set them directly on the plotting area on the web page
        #     # hovermode = closest shows the values for the hover point, otherwise by default ('compare' mode) you only see one coordinate
        #     combined_figure.layout.update(dict(title='Combined scores', autosize=False, width=600, height=plotheight, hovermode='closest', showlegend=False))
        #     py.plot(combined_figure, filename=os.path.join(self.plots_folder, 'heatmap_combined_data.html'), auto_open=False, show_link=False, include_plotlyjs=True)
