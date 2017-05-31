def includeme(config):
    config.add_static_view('static', 'static', cache_max_age=3600)
    config.add_route('home', '/')
    config.add_route('projects', '/projects')
    config.add_route('project_add', '/projects/add')
    config.add_route('project_view', '/projects/{projectid}')
    config.add_route('project_edit', '/projects/{projectid}/edit')
    config.add_static_view('deform_static', 'deform:static/')
