def includeme(config):
    config.add_static_view('static', 'static', cache_max_age=3600)
    config.add_static_view('generated', 'generated', cache_max_age=3600)
    config.add_static_view('deform_static', 'deform:static/')

    # home page: projects table and upload project layout file
    config.add_route('home', '/')

    # individual project page: view (overvire, results & samples), add comments and load data
    config.add_route('project_view', '/project/{projectid}')
    config.add_route('project_edit', '/project/{projectid}/edit')
    config.add_route('project_sequence', '/project/{projectid}/sequence')
