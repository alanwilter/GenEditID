import os
import uuid
import shutil

import colander
import deform.widget
import datetime

from pyramid.httpexceptions import HTTPFound
from pyramid.view import view_config

from dnascissors.model import Project


# See http://docs.pylonsproject.org/projects/pyramid/en/latest/quick_tutorial/forms.html

# File uploads: http://docs.pylonsproject.org/projects/pyramid-cookbook/en/latest/forms/file_uploads.html

class LoaderViews(object):
    
    def __init__(self, request):
        self.request = request
        self.dbsession = request.dbsession
    
    @view_config(route_name="load_layout", renderer="../templates/loader/layoutsetup.pt")
    def load_layout(self):
        
        if 'submit' in self.request.params:
            
            filename = self.request.POST['layoutfile'].filename
            filedata = self.request.POST['layoutfile'].file
            
            print("Uploaded = %s" % filename)
            
            file_path = os.path.join('uploads/', '%s.xlsx' % uuid.uuid4())
    
            temp_file_path = file_path + '~'
    
            try:
                filedata.seek(0)
                with open(temp_file_path, 'wb') as output_file:
                    shutil.copyfileobj(filedata, output_file)
                
                os.rename(temp_file_path, file_path)
                
                # Now do the load, after cleaning project.
                
                statinfo = os.stat(file_path)
                
                print("Uploaded a file of {:d} bytes".format(statinfo.st_size))
                
            except Exception as e:
                return dict(title="Load Layout File")
            finally:
                try:
                    os.remove(temp_file_path)
                except OSError:
                    # Do nothing.
                    pass
                    
                try:
                    os.remove(file_path)
                except OSError:
                    # Do nothing.
                    pass

            url = self.request.route_url('projects')
            return HTTPFound(url)
        
        return dict(title="Load Layout File")
