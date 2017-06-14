import colander
import deform.widget
import datetime

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from dnascissors.model import Project, ExperimentLayout

from webapp.plots.plotter import Plotter

class ExperimentViews(object):
    
    def __init__(self, request):
        self.request = request
        self.dbsession = request.dbsession
    
    @view_config(route_name="experiment_view", renderer="../templates/experiment/viewlayout.pt")
    def viewLayout(self):
        layoutid = self.request.matchdict['layoutid']
        
        layout = self.dbsession.query(ExperimentLayout).filter(ExperimentLayout.id == layoutid).one()
        
        plotter = Plotter()
        
        return dict(layout=layout,
                    title="Gene Editing Experiment %s" % layout.geid,
                    cellgrowthplot=plotter.growth_plot(self.dbsession, layout.project.geid, layout.geid),
                    proteinabundanceplot=plotter.abundance_plot(self.dbsession, layout.project.geid, layout.geid))
        
