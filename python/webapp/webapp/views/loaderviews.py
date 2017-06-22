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
from dnascissors.loader import LayoutLoader, ExistingEntityException


# See http://docs.pylonsproject.org/projects/pyramid/en/latest/quick_tutorial/forms.html

# File uploads: http://docs.pylonsproject.org/projects/pyramid-cookbook/en/latest/forms/file_uploads.html

class LoaderViews(object):
    
    def __init__(self, request):
        self.logger = logging.getLogger('webapp')
        self.request = request
        self.dbsession = request.dbsession
    
    @view_config(route_name="load_layout", renderer="../templates/loader/layoutsetup.pt")
    def load_layout(self):
        
        return_map = dict(title="Load Layout File")
        
        if 'submit' in self.request.params:
            
            filename = self.request.POST['layoutfile'].filename
            filedata = self.request.POST['layoutfile'].file
            
            clean_existing = False
            try:
                clean_existing = self.request.POST['blat']
            except KeyError:
                pass
            
            self.logger.info("Uploaded = %s" % filename)
            
            file_path = os.path.join('uploads/', '%s.xlsx' % uuid.uuid4())
    
            temp_file_path = file_path + '~'
    
            try:
                filedata.seek(0)
                with open(temp_file_path, 'wb') as output_file:
                    shutil.copyfileobj(filedata, output_file)
                
                os.rename(temp_file_path, file_path)
                
                statinfo = os.stat(file_path)
                
                self.logger.info("Uploaded a file of {:d} bytes".format(statinfo.st_size))
                
                # Now do the load, after cleaning project.
                
                loader = LayoutLoader(self.dbsession, file_path)
                
                try:
                    loader.load_all(clean_existing)

                    url = self.request.route_url('projects')
                    return HTTPFound(url)
                
                except ExistingEntityException as e:
                    return_map['error'] = e.message
                
            except Exception as e:
                self.logger.error("Have an unexpected error while creating project: {}".format(e))
                return_map['error'] = str(e)
            
            finally:
                try:
                    os.remove(temp_file_path)
                except OSError:
                    pass
                    
                try:
                    os.remove(file_path)
                except OSError:
                    pass
        
        return return_map
