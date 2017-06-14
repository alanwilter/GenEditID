import re
import colorlover

import plotly
import plotly.graph_objs as go
import plotly.offline as py

import sqlalchemy
import logging

from datetime import datetime
from pandas.core.frame import DataFrame

from dnascissors.config import cfg
from dnascissors.model import *

from naturalsort import natural_sort


class Plotter:
    
    def create_classifier(self, cell_line, clone, guide, well_content):
        
        parts = []
        
        if cell_line:
            parts.append(cell_line.name)
        
        if clone:
            parts.append(clone.name)
        
        if guide:
            parts.append(guide.name)
        
        if well_content:
            content = well_content.content_type
            if content and content not in ['empty', 'sample']:
                parts.append(content)
        
        return " ".join(parts)
        

    def growth_plot(self, dbsession, projectid, plateid):
        
        query = dbsession.query(CellGrowth).\
                join(CellGrowth.well).\
                join(Well.well_content).\
                join(Well.experiment_layout).\
                join(ExperimentLayout.project).\
                outerjoin(WellContent.clone).\
                join(Clone.cell_line).\
                filter(Project.id == projectid).\
                filter(ExperimentLayout.geid == plateid)
                #filter(CellGrowth.confluence_percentage != None)

        all_growth = query.all()
        
        # Assemble list of lists to create a data frame.
        
        data = []
        
        for growth in all_growth:
            well = growth.well
            clone = well.well_content.clone
            cell_line = clone.cell_line
            layout = well.experiment_layout
            
            #guideNames = [lambda g: g.name for g in well.well_content.guides]
            #guides = ', '.join(guideNames)
            guides = ', '.join([g.name for g in well.well_content.guides])
            
            position = "%s%d" % (well.row, well.column)
            
            for guide in well.well_content.guides:
                target = guide.target
            
                row = [ self.create_classifier(cell_line, clone, guide, well.well_content),
                        layout.geid, target.name, guide.name, position, growth.hours, growth.confluence_percentage ]
            
                data.append(row)
        
        # Pandas data frame
        
        data = DataFrame(data, columns=['content', 'plate', 'target', 'guide', 'well', 'elapsed', 'confluence'])
        
        sorted_wells = data.well.unique().tolist()
        natural_sort(sorted_wells)
        
        #unique_content = data.content.unique()
        
        #colours = colorlover.scales[str(len(sorted_wells))]['div']
        #colours = colorlover.scales['3']['div']['RdYlBu']
        #print(colours)
        #colours = colorlover.interp(colours, len(unique_content))
        #print(colours)
        
        #colour_map = dict()
        #for i in range(0, len(unique_content)):
        #    colour_map[unique_content[i]] = colours[i]
        
        # Need to assemble several plot objects for each line.
        
        plots = []
        
        for well in sorted_wells:
        
            welldata = data[data.well == well]
            
            # See https://stackoverflow.com/a/16729808
            content = welldata.iloc[0]['content']
            
            #colour = colour_map[content]
            colour = 'blue'
        
            plots.append(
                go.Scatter(
                    mode='lines',
                    line=dict(color=colour),
                    x=welldata.elapsed,
                    y=welldata.confluence,
                    name="%s (%s)" % (well, content),
                    hoverinfo='none'
                )
            )
        
        layout = go.Layout(
            title="Growth of Plate %s" % plateid,
            xaxis=dict(title="Time (h)"),
            yaxis=dict(title="Confluence (%)", range=[0, 100])
        )

        figure = go.Figure(data=plots, layout=layout)
        
        py.plot(figure, filename="growth_plate_%s.html" % plateid)


    def abundance_plot(self, dbsession, projectid, plateid):
        
        query = dbsession.query(ProteinAbundance).\
                join(ProteinAbundance.well).\
                join(Well.well_content).\
                join(Well.experiment_layout).\
                join(ExperimentLayout.project).\
                outerjoin(WellContent.clone).\
                join(Clone.cell_line).\
                filter(Project.id == projectid).\
                filter(ExperimentLayout.geid == plateid).\
                filter(WellContent.content_type != 'background')

        all_abundance = query.all()
        
        # Assemble list of lists to create a data frame.
        
        data = []
        
        for abundance in all_abundance:
            well = abundance.well
            clone = well.well_content.clone
            cell_line = clone.cell_line
            layout = well.experiment_layout
            
            position = "%s%d" % (well.row, well.column)
            
            for guide in well.well_content.guides:
                target = guide.target
            
                row = [ self.create_classifier(cell_line, clone, guide, well.well_content),
                        layout.geid, position, abundance.intensity_channel_800 / abundance.intensity_channel_700 ]
            
                data.append(row)
        
        # Pandas data frame
        
        data = DataFrame(data, columns=['content', 'plate', 'well', 'channelratio'])

        unique_content = data.content.unique().tolist()
        natural_sort(unique_content)
        
        #colours = colorlover.scales['3']['div']['RdYlBu']
        #print(colours)
        #print(unique_content)
        #if len(unique_content) > 3:
        #    colours = colorlover.interp(colours, len(unique_content) + 1)
        #print(colours)
        
        #colour_map = dict()
        #for i in range(0, len(unique_content)):
        #    colour_map[unique_content[i]] = colours[i]
        
        # Need to assemble several plot objects for each line.
        
        colour = 'blue'
    
        plots = []
        
        for content in unique_content:
        
            contentdata = data[data.content == content]
            
            plots.append(
                go.Scatter(
                    mode='markers',
                    line=dict(color=colour),
                    x=contentdata.content,
                    y=contentdata.channelratio,
                    name=content,
                    hoverinfo='none'
                )
            )

        layout = go.Layout(
            title="Protein Abundance of Plate %s" % plateid,
            xaxis=dict(title="Cell Line", ticks=False, fixedrange=True),
            yaxis=dict(title="Relative protein abundance")
        )

        figure = go.Figure(data=plots, layout=layout)
        
        py.plot(figure, filename="abundance_plate_%s.html" % plateid)

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
        
        plotter.growth_plot(session, 1, 'GEP00001_01')
        
        plotter.abundance_plot(session, 1, 'GEP00001_02')
    
    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
