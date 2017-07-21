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

    def __init__(self):
        self.include_js = False

    def calculate_percentage_plots(self, dfgroup, grouping_variable, number_loop = -1, variant_caller = None):
        plots =[]
        marker_symbol = ['circle', 'triangle-up', 'cross', 'hash']
        nloop = number_loop + 1
        msym = 'circle' if variant_caller == None else marker_symbol[nloop]
        for guide_name, grouped_data in dfgroup:
            grouped_data_byvar = grouped_data.groupby([grouping_variable]).size()
            grouped_data_byvar_percent = grouped_data_byvar*100 / grouped_data_byvar.sum()
            htext = []
            for length_gdbp in range(len(grouped_data_byvar_percent)):
                if variant_caller != None:
                    hovertext = [variant_caller] #add other elements to this list to display when hovering
                else:
                    hovertext = []
                    hovertext = ' '.join(hovertext)
                    htext.append(hovertext)
            plots.append(
                go.Scatter(
                     x = grouped_data_byvar_percent.index.tolist(),
                     y = grouped_data_byvar_percent.tolist(),
                     name = guide_name,
                     mode = 'markers',
                     marker = {'symbol': msym},
                     text = htext
                )
            )
        return(plots)

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

    def growth_plot(self, dbsession, projectid, file=None):
        # See https://stackoverflow.com/questions/21114830/query-to-check-if-size-of-collection-is-0-or-empty-in-sqlalchemy
        # Get all wells for the project where there is growth information.
        query = dbsession.query(Well)\
                         .join(Well.well_content)\
                         .join(Well.experiment_layout)\
                         .join(ExperimentLayout.project)\
                         .filter(Project.geid == projectid)\
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
        if not file:
            output_type = "div"
            file = "cell_growth_{}.html".format(projectid)
        return py.plot(figure, filename=file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)

    def abundance_plot(self, dbsession, projectid, file=None):
        # Get all wells for the project where there is protein abundance information.
        query = dbsession.query(Well)\
                         .join(Well.well_content)\
                         .join(Well.experiment_layout)\
                         .join(ExperimentLayout.project)\
                         .filter(Project.geid == projectid)\
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
        if not file:
            output_type = "div"
            file = "protein_abundance_{}.html".format(projectid)
        return py.plot(figure, filename=file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)

    def zygosity_plot(self, dbsession, projectid, file=None):
        # we need to filter by gDNA, because there is a sample in project1 (GE-P6B4-G)
        # that is gDNA only (it was sent for sequencing only as gDNA, without a 'fixed cells' counterpart)
        query = dbsession.query(Well)\
                .join(Well.sequencing_library_contents)\
                .join(Well.experiment_layout)\
                .join(ExperimentLayout.project)\
                .filter(Project.geid == projectid)\
                .filter(SequencingLibraryContent.dna_source != 'gDNA')

        wells = query.all()
        if len(wells) == 0:
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
            # print(well.sequencing_library_contents[0].sequencing_sample_name, guide_name, mutation_zygosity)
            zygosities.append(mutation_zygosity)
            guides.append(guide_name)
            
        # convert 'results' to pandas dataframe and group by 'guides'
        df = pandas.DataFrame({'guides': guides, 'zygosities': zygosities})
        dfgroup = df.groupby(['guides'])
        
        plots = self.calculate_percentage_plots(dfgroup, 'zygosities')
        # order x axis values
        categories = ['wt', 'homo', 'smut', 'dmut', 'iffy']
        #print(categories)
        layout = go.Layout(
            title='Zygosities',
            xaxis={'categoryorder': 'array', 'categoryarray': categories},
            yaxis={'title': '% of submitted samples per guide'}
        )
        # plot                       
        figure = go.Figure(data=plots, layout=layout)
        output_type = "file"
        if not file:
            output_type = "div"
            file = "plot_zygosity_{}.html".format(projectid)
        return py.plot(figure, filename=file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)


    def typeofvariant_plot(self, dbsession, projectid, file_INDELs=None, file_SNVs=None):
        
        query = dbsession.query(VariantResult)\
                         .join(VariantResult.sequencing_library_content)\
                         .join(SequencingLibraryContent.well)\
                         .join(Well.well_content)\
                         .join(Well.experiment_layout)\
                         .join(ExperimentLayout.project)\
                         .filter(Project.geid == projectid)\
                         .filter(SequencingLibraryContent.dna_source != 'gDNA')\
                         .filter(VariantResult.allele_fraction > 0.1)
        wells = query.all()               
        if len(wells) == 0:
            return None
        
        guide = []
        typemutation = []
        typevariant = []
        typecaller = []
        for i in wells:
            well = i.sequencing_library_content.well
            guide_name = 'none'
            
            if well.well_content.guides:
                guide_name = well.well_content.guides[0].name
           # print(guide_name, i.consequence, i.variant_type, i.variant_caller)
            
            guide.append(guide_name)
            typemutation.append(i.consequence)
            typevariant.append(i.variant_type)
            typecaller.append(i.variant_caller)
        
        #print(len(wells))


        #-------- calculate percentages per guide in a pandas dataframe
        # convert 'results' to pandas dataframe and group by 'variants'
        df = pandas.DataFrame({'variants': typevariant, 'caller': typecaller, 'guides': guide, 'mutation': typemutation})
        dfgroup_typevariant = df.groupby(['variants'])
        # create independent datasets (one for indels, one for snvs) for plotting
        dfgroup_typevariant_INDEL = dfgroup_typevariant.get_group('INDEL')  #[i for i in dfgroup_typevariant][0][1]
        dfgroup_typevariant_SNV = dfgroup_typevariant.get_group('SNV')      #[i for i in dfgroup_typevariant][1][1]
        
        # calculate percentages of mutation types and create bar plot 'data' dictionary
        # This will produce a plot with grouped legends. Annoyingly there no feature at date 20170710 to add 
        #  titles to the legend groups.
        #  It's been suggested to use layout annotations as a workaround: https://github.com/plotly/plotly.js/issues/689
        #  I assume that the top legend group is the first variant caller.
        plot_indels = []
        nloop = -2
        for caller, gdata in dfgroup_typevariant_INDEL.groupby(['caller']):
            nloop += 1
            plotindels = self.calculate_percentage_plots(gdata.groupby(['guides']), 'mutation', number_loop = nloop, variant_caller = caller)
            plot_indels.extend(plotindels)
            
        plots = []
        nloop = -2
        for caller, gdata in dfgroup_typevariant_SNV.groupby(['caller']):
            nloop += 1
            plots_SNVs = self.calculate_percentage_plots(gdata.groupby(['guides']), 'mutation', number_loop = nloop, variant_caller = caller)
         
        # layouts
        layout_INDELs = go.Layout(
            title = 'Type of mutation (INDELS)',
            yaxis = {'title': '% of submitted samples per guide'}    
        )
        
        layout_SNVs = go.Layout(
            title = 'Type of mutation (SNVs)',
            yaxis = {'title': '% of submitted samples per guide'}    
        )
        
        # plot
        figure_INDELs = go.Figure(data=plot_indels, layout=layout_INDELs)
        #figure_SNVs = go.Figure(data=plots_SNVs, layout=layout_SNVs)
        output_type = "file"
        
        if not file_INDELs:
            output_type = "div"
            file_INDELs = "plotNGS_INDELs_{}.html".format(projectid)
            
        return py.plot(figure_INDELs, filename=file_INDELs, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)
                      
        if not file_SNVs:
            output_type = "div"
            file_INDELs = "plotNGS_SNVs_{}.html".format(projectid)
            
        #return py.plot(figure_SNVs, filename=file_SNVs, auto_open=False, show_link=False,
         #              include_plotlyjs=self.include_js, output_type=output_type)
       # return py.plot(dict(data = data_INDELs, layout = layout_INDELs), filename = 'plotNGS_indels.html')             
            # not sure why, but when the script is sourced, both plots are the SNV plot
        
        
        # the ideal situation would be to plot both datasets side by side
         # pyoff.plot(dict(data = data_INDELs + data_SNVs, layout = layout_INDELs)) #say we use a different layout for each dataset
         # however the guides appear duplicated. As far as I have seen, you can't have 
         # the same legent for two different plots in plotly (ironically you can have it in if you use ggplotly(R:ggplot2)... )
                
    def indellengths_plot(self, dbsession, projectid, file = None):
        
        query = dbsession.query(VariantResult)\
             .join(VariantResult.sequencing_library_content)\
             .join(SequencingLibraryContent.well)\
             .join(Well.experiment_layout)\
             .join(ExperimentLayout.project)\
             .filter(Project.geid == projectid)\
             .filter(SequencingLibraryContent.dna_source != 'gDNA')\
             .filter(VariantResult.allele_fraction > 0.1)
        wells = query.all()               
        if len(wells) == 0:
            return None
        guides = []
        allindellengths = []
        typecallers = []
        
        for i in wells:
            well = i.sequencing_library_content.well
            layout = well.experiment_layout
            guide_name = 'none'
            
            if well.well_content.guides: 
                guide_name = well.well_content.guides[0].name
            # print(i.sequencing_library_content.sequencing_sample_name, guide_name, i.indel_length, i.variant_caller)
            guides.append(guide_name)    
            allindellengths.append(i.indel_length)
            typecallers.append(i.variant_caller)
        
        
        #-------- calculate percentages per guide in a pandas dataframe
        # convert 'results' to pandas dataframe and group by 'guides'
        df = pandas.DataFrame({'caller': typecallers, 'guides': guides, 'indellengths': allindellengths})
        # calculate percentages of indel lengths and create bar plot 'data' dictionary
        plots = []
        nloop = -2
        for caller, gdata in df.groupby(['caller']):
            nloop += 1
            plot_indellenghts = self.calculate_percentage_plots(gdata.groupby(['guides']), 'indellengths', number_loop = nloop, variant_caller = caller)
            plots.extend(plot_indellenghts)
        
        
        print(plots)
        
        
        layout = go.Layout(
            title = 'Indel lengths',
            yaxis = {'title': '% of submitted samples per guide'}    
        )
        
        # plot
        figure = go.Figure(data=plots, layout=layout)
        #figure_SNVs = go.Figure(data=plots_SNVs, layout=layout_SNVs)
        output_type = "file"
        
        if not file:
            output_type = "div"
            file = "plotNGS_indellengths_{}.html".format(projectid)
            
        return py.plot(figure, filename=file, auto_open=False, show_link=False,
                       include_plotlyjs=self.include_js, output_type=output_type)
                      
def main():
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M')
    engine = sqlalchemy.create_engine(cfg['DATABASE_URI'])
    Base.metadata.bind = engine
    DBSession = sqlalchemy.orm.sessionmaker(bind=engine)
    session = DBSession()
    try:
        plotter = Plotter()
        plotter.include_js = True
        plotter.growth_plot(session, 'GEP00001', "growth.html")
        plotter.abundance_plot(session, 'GEP00001', "abundance.html")
        plotter.zygosity_plot(session, 'GEP00001', "plotNGS_zygosity.html")
        plotter.typeofvariant_plot(session, 'GEP00001', "plotNGS_INDELs.html")
        plotter.indellengths_plot(session, 'GEP00001', "plotNGS_indellengths.html" )
  #      plotter.typeofvariant_plot(session, 'GEP00001', "plotNGS_SNVs.html")        
    
    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
