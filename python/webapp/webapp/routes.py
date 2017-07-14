def includeme(config):
    config.add_static_view('static', 'static', cache_max_age=3600)
    config.add_static_view('generated', 'generated', cache_max_age=3600)

    config.add_route('home', '/')

    config.add_route('load_layout', '/loadlayout')

    config.add_route('projects', '/project')
    config.add_route('project_view', '/project/{projectid}')
    config.add_route('project_edit', '/project/{projectid}/edit')
    config.add_route('project_plot', '/project/{projectid}/plot')

    config.add_route('experiment_view', '/experiment/{layoutid}')

    config.add_route('plate_load', '/plate/{plateid}')
    config.add_route('plate_icw', '/plate/{plateid}/icw')
    config.add_route('plate_incu', '/plate/{plateid}/incucyte')

    config.add_static_view('deform_static', 'deform:static/')
