import re
import colorlover

import plotly
#import plotly.plotly as py
import plotly.graph_objs as go
import plotly.offline as py

import sqlalchemy
import logging

from datetime import datetime
from pandas.core.frame import DataFrame

from dnascissors.config import cfg
from dnascissors.model import *

from naturalsort import natural_sort

#from plotly.offline import download_plotlyjs, init_notebook_mode, plot, iplot

class Plotter:
    
    def create_classifier(self, cell_line, clone, guide, content):
        
        parts = []
        
        if cell_line:
            parts.append(cell_line.name)
        
        if clone:
            parts.append(clone.name)
        
        if guide:
            parts.append(guide.name)
        
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
            
                row = [ self.create_classifier(cell_line, clone, guide, well.well_content.content_type),
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
        
        plotter.growth_plot(session, 1, 'GEP00001_01_incu')
    
    except Exception as e:
        logging.exception(e)
    finally:
        session.close()

if __name__ == '__main__':
    main()
