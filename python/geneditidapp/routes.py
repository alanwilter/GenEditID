def includeme(config):
    config.add_static_view('static', 'static', cache_max_age=3600)
    config.add_static_view('generated', 'generated', cache_max_age=3600)
    config.add_static_view('deform_static', 'deform:static/')

    # home page: create new project and projects table
    config.add_route('home', '/')

    # help page: getting started
    config.add_route('help', '/help')

    # individual project page: configure & view results
    config.add_route('project', '/project/{gepid}')
    #config.add_route('project_edit', '/project/{projectid}/edit')
