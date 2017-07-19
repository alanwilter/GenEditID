import logging
import os
import shutil
import uuid

from pyramid.view import view_config

from dnascissors.model import Plate, ProteinAbundance, CellGrowth
from dnascissors.loader import ProteinAbundanceLoader, CellGrowthLoader


class ExperimentViews(object):

    def __init__(self, request):
        self.logger = logging.getLogger('webapp')
        self.request = request
        self.dbsession = request.dbsession

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
            if self.request.POST["icwfile"] is not None:
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
            if self.request.POST["incufile"] is not None:
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
