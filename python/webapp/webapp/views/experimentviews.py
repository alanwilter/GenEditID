import colander
import deform.widget
import datetime
import logging
import os
import shutil
import uuid

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from dnascissors.model import Project, ExperimentLayout, Plate, ProteinAbundance, CellGrowth
from dnascissors.loader import ProteinAbundanceLoader, CellGrowthLoader, ExistingEntityException

from webapp.plots.plotter import Plotter

class ExperimentViews(object):
    
    def __init__(self, request):
        self.logger = logging.getLogger('webapp')
        self.request = request
        self.dbsession = request.dbsession
    
    @view_config(route_name="experiment_view", renderer="../templates/experiment/viewlayout.pt")
    def view_layout(self):
        layoutid = self.request.matchdict['layoutid']
        
        layout = self.dbsession.query(ExperimentLayout).filter(ExperimentLayout.id == layoutid).one()
        
        plotter = Plotter()
        
        headers = [ "Plate", "Well", "Sample", "Barcode", "Score", "Protein abundance to WT", "Protein abundance to KO",\
                    "Growth slope to WT", "Growth slope to KO", "Variant Type", "Symbol", "Confidence", "Allele fraction", "Alleles" ]
        rows = []
        
        for well in layout.wells:
            
            guide = None
            
            anyDone = False
            
            for slc in well.sequencing_library_contents:
                for vc in slc.variant_results:
                    
                    row = []
                    row.append(layout.geid)
                    row.append("{:s}{:02}".format(well.row, well.column))
                    row.append(slc.sequencing_sample_name)
                    row.append(slc.sequencing_barcode)
                    row.append("-")
                    row.append("-")
                    row.append("-")
                    row.append("-")
                    row.append("-")
                    row.append(vc.consequence)
                    row.append(vc.gene)
                    row.append("-")
                    row.append("{0:.3f}".format(vc.allele_fraction))
                    row.append(vc.alleles)
                    
                    rows.append(row)
                    
                    anyDone = True
            
            """
            if not anyDone:
                    row = [None for i in range(len(headers))]
                    row[0] = "{:s}{:02}".format(well.row, well.column)
                    row[5] = guide
                    rows.append(row)
            """
            
        return dict(layout=layout,
                    title="Gene Editing Experiment %s" % layout.geid,
                    column_headers=headers,
                    rows=rows,
                    cellgrowthplot=plotter.growth_plot(self.dbsession, layout.project.geid, layout.geid),
                    proteinabundanceplot=plotter.abundance_plot(self.dbsession, layout.project.geid, layout.geid))
    
    @view_config(route_name="plate_load", renderer="../templates/experiment/loadplate.pt")
    def load_plate(self):
        plateid = self.request.matchdict['plateid']
        
        plate = self.dbsession.query(Plate).filter(Plate.id == plateid).one()
        
        return dict(plate=plate,
                    title="Gene Editing Plate %s" % plate.geid,
                    icwplate=self._has_icw(plate.id),
                    incuplate=self._has_incu(plate.id))
    
    @view_config(route_name="plate_icw", renderer="../templates/experiment/loadplate.pt")
    def load_icw(self):
        plateid = self.request.matchdict['plateid']

        plate = self.dbsession.query(Plate).filter(Plate.id == plateid).one()

        return_map = dict(plate=plate,
                          title="Gene Editing Plate %s" % plate.geid)
        
        if 'submiticw' in self.request.params:
            
            if self.request.POST["icwfile"] != None:
                
                file_path = None
                
                try:
                    file_path = self._upload("icwfile", ".csv")
                    
                    if file_path:
                        
                        loader = ProteinAbundanceLoader(self.dbsession, file_path, plate.geid)
                        
                        self._remove_existing_rows(plateid)
                        
                        loader.load()
                    
                except Exception as e:
                    self.logger.error("Have an unexpected error while uploading ICW: {}".format(e))
                    raise
                
                finally:
                    if file_path:
                        try:
                            os.remove(file_path)
                        except OSError:
                            pass
        
        return_map['icwplate'] = self._has_icw(plate.id)
        return_map['incuplate'] = self._has_incu(plate.id)
        
        return return_map

    @view_config(route_name="plate_incu", renderer="../templates/experiment/loadplate.pt")
    def load_incucyte(self):
        plateid = self.request.matchdict['plateid']

        plate = self.dbsession.query(Plate).filter(Plate.id == plateid).one()
        
        return_map = dict(plate=plate,
                          title="Gene Editing Plate %s" % plate.geid)
        
        if 'submitincu' in self.request.params:
            
            if self.request.POST["incufile"] != None:
                
                file_path = None
                
                try:
                    file_path = self._upload("incufile", ".txt")
                    
                    if file_path:
                        
                        loader = CellGrowthLoader(self.dbsession, file_path, plate.geid)
                        
                        self._remove_existing_rows(plateid)
                        
                        loader.load()
                    
                except Exception as e:
                    self.logger.error("Have an unexpected error while uploading Incucyte: {}".format(e))
                    raise
                
                finally:
                    if file_path:
                        try:
                            os.remove(file_path)
                        except OSError:
                            pass

        return_map['icwplate'] = self._has_icw(plate.id)
        return_map['incuplate'] = self._has_incu(plate.id)
        
        return return_map

    def _upload(self, property, suffix): 
        
        filename = self.request.POST[property].filename
        filedata = self.request.POST[property].file
        
        if not filedata:
            return None
            
        self.logger.debug("Uploaded = %s" % filename)
            
        file_path = os.path.join('uploads/', "{}{}".format(uuid.uuid4(), suffix))

        temp_file_path = file_path + '~'

        try:
            filedata.seek(0)
            with open(temp_file_path, 'wb') as output_file:
                shutil.copyfileobj(filedata, output_file)
            
            os.rename(temp_file_path, file_path)
            
            statinfo = os.stat(file_path)
            
            self.logger.info("Uploaded a file of {:d} bytes".format(statinfo.st_size))
            
        finally:
            try:
                os.remove(temp_file_path)
            except OSError:
                pass
        
        return file_path

    def _has_icw(self, plateid):
        
        rows = self.dbsession.query(ProteinAbundance).filter(ProteinAbundance.plate_id == plateid).count()
        
        return rows > 0
    
    def _has_incu(self, plateid):
        
        rows = self.dbsession.query(CellGrowth).filter(CellGrowth.plate_id == plateid).count()
        
        return rows > 0

    def _remove_existing_rows(self, plateid):
        
        # Removes both cell growth and protein abundance since a plate can't be for both.
        
        self.dbsession.query(ProteinAbundance).filter(ProteinAbundance.plate_id == plateid).delete()
        self.dbsession.query(CellGrowth).filter(CellGrowth.plate_id == plateid).delete()
        self.dbsession.flush()
        