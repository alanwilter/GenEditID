import os
import uuid
import shutil
import logging

import colander
import deform.widget
import datetime

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from dnascissors.model import Project
from dnascissors.loader import LayoutLoader, ProteinAbundanceLoader, ExistingEntityException


# See http://docs.pylonsproject.org/projects/pyramid/en/latest/quick_tutorial/forms.html

# File uploads: http://docs.pylonsproject.org/projects/pyramid-cookbook/en/latest/forms/file_uploads.html

class LoaderViews(object):
    
    def __init__(self, request):
        self.logger = logging.getLogger('webapp')
        self.request = request
        self.dbsession = request.dbsession
    
    @view_config(route_name="load_layout", renderer="../templates/loader/layoutsetup.pt")
    def load_layout(self):
        
        return_map = dict(title="Load Layout File", clash=False)
        
        if 'submit' in self.request.params:
            
            clean_existing = False
            try:
                clean_existing = self.request.POST['blat']
            except KeyError:
                pass
            
            file_path = None
            
            try:
                file_path = self._upload("layoutfile", ".xlsx")
                
                # Now do the load, after cleaning project.
                
                loader = LayoutLoader(self.dbsession, file_path)
                
                try:
                    loader.load_all(clean_existing)

                    url = self.request.route_url('projects')
                    return HTTPFound(url)
                
                except ExistingEntityException as e:
                    return_map['error'] = e.message
                    return_map['clash'] = e
                    return_map['layoutfile'] = self.request.POST['layoutfile']
                
            except Exception as e:
                self.logger.error("Have an unexpected error while creating project: {}".format(e))
                return_map['error'] = str(e)
            
            finally:
                if file_path:
                    try:
                        os.remove(file_path)
                    except OSError:
                        pass
        
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
