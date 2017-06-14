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
        
        # SQLAlchemy query. Roughly equates to:
        """
            select distinct
                cl.name as cellline,
                c.name as clone,
                g.name as guide,
                wc.content_type as type,
                p.geid as plate,
                concat(w.row, w.column) as well,
                t.name as target,
                gr.hours as elapsed,
                gr.confluence_percentage as confluence
            from well w
                inner join growth gr on gr.well_id = w.id
                inner join well_content wc on w.well_content_id = wc.id
                inner join experiment_layout el on w.experiment_layout_id = el.id
                inner join plate p on p.experiment_layout_id = el.id
                left join clone c on wc.clone_id = c.id
                inner join cell_line cl on c.cell_line_id = c.id
                left join guide_well_content_association gwca on gwca.well_content_id = wc.id
                inner join guide g on gwca.guide_id = g.id
                inner join target t on g.target_id = t.id
                inner join project on t.project_id = project.id
            where
                gr.confluence_percentage is not null
                and project.id = :projectid
                and plate.id = :plateid
        """
        

        
        query = dbsession.query(CellGrowth).\
                join(CellGrowth.plate).\
                join(CellGrowth.well).\
                join(Well.well_content).\
                join(Well.experiment_layout).\
                outerjoin(WellContent.clone).\
                join(Clone.cell_line).\
                outerjoin(WellContent.guides).\
                join(Guide.target).\
                join(Target.project).\
                filter(Project.id == projectid).\
                filter(Plate.geid == plateid).\
                filter(CellGrowth.confluence_percentage != None)

        all_growth = query.all()
        
        # Assemble list of lists to create a data frame.
        
        data = []
        
        for growth in all_growth:
            well = growth.well
            clone = well.well_content.clone
            cell_line = clone.cell_line
            plate = growth.plate
            
            position = "%s%d" % (well.row, well.column)
            
            for guide in well.well_content.guides:
                target = guide.target
            
                row = [ self.create_classifier(cell_line, clone, guide, well.well_content),
                        plate.geid, target.name, guide.name, position, growth.hours, growth.confluence_percentage ]
            
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
            title="Growth of Plate %s" % plate.geid,
            xaxis=dict(title="Time (h)"),
            yaxis=dict(title="Confluence (%)", range=[0, 100])
        )

        figure = go.Figure(data=plots, layout=layout)
        
        py.plot(figure, filename="growth_plate_%s.html" % plateid)


    def abundance_plot(self, dbsession, projectid, plateid):
        
        # The dplyr R version:
        """
            left_join(well, abundance, by="well_id") %>%
            left_join(., well_content, by="well_content_id") %>%
            left_join(., clone, by="clone_id") %>%
            left_join(., cell_line, by="cell_line_id") %>%
            left_join(., plate, by="plate_id") %>%
            left_join(., guide_well_content_association, by="well_content_id") %>%
            left_join(., guide, by="guide_id") %>%
            left_join(., target, by="target_id") %>%
            collect %>%
            filter(!is.na(content_type)) %>%
            mutate(Well = paste0(row, column)) %>%
            mutate(classifier=create_classifier(cell_line_name, clone_name, guide_name, content_type)) %>%
            select(classifier, plate_id, Well, intensity_channel_700, intensity_channel_800) %>%
            rename(Content=classifier, Plate = plate_id, Channel700=intensity_channel_700, Channel800=intensity_channel_800) %>%
            filter(!grepl('background', Content)) %>%            #the background control (no antibodies is not meaningful for proteini amount)
            mutate('ratio800to700' = Channel800/Channel700) %>%  #calculate relative protein abundance, ratio 800 to 700
            mutate(Plate = as.factor(Plate), Well = as.factor(Well))  #change Plate and position to factors, so grouping in ggplot can work
        """
        
        query = dbsession.query(ProteinAbundance).\
                join(ProteinAbundance.plate).\
                join(ProteinAbundance.well).\
                join(Well.well_content).\
                join(Well.experiment_layout).\
                outerjoin(WellContent.clone).\
                join(Clone.cell_line).\
                outerjoin(WellContent.guides).\
                join(Guide.target).\
                join(Target.project).\
                filter(Project.id == projectid).\
                filter(Plate.geid == plateid).\
                filter(WellContent.content_type != 'background')

        all_abundance = query.all()
        
        # Assemble list of lists to create a data frame.
        
        data = []
        
        for abundance in all_abundance:
            well = abundance.well
            clone = well.well_content.clone
            cell_line = clone.cell_line
            plate = abundance.plate
            
            position = "%s%d" % (well.row, well.column)
            
            for guide in well.well_content.guides:
                target = guide.target
            
                row = [ self.create_classifier(cell_line, clone, guide, well.well_content),
                        plate.geid, position, abundance.intensity_channel_800 / abundance.intensity_channel_700 ]
            
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
            
            print(content)
            print(contentdata)
            
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
            title="Protein Abundance of Plate %s" % plate.geid,
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
        
        #plotter.growth_plot(session, 1, 'GEP00001_01_incu')
        
        plotter.abundance_plot(session, 1, 'GEP00001_02_ICW')
    
    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
