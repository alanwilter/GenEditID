import logging

from pyramid.response import Response
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPFound


class HelpViews(object):
    def __init__(self, request):
        self.logger = logging.getLogger(__name__)

    @view_config(route_name='help', renderer='../templates/help.pt')
    def ge_help_page(self):
        return_map = {'title': 'GenEditID',
                      'error': False,
                      'info': False,
                     }
        return return_map
