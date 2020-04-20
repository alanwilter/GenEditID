from setuptools import setup, find_packages

version = '0.1'

print("""
-----------------------------------
 Installing GeneditID version {}
-----------------------------------
""".format(version))

install_requires = [
    # configuration
    'pyyaml',
    # database
    # 'psycopg2',  # for production only
    'sqlalchemy',
    # data
    'xlrd',
    'pandas',
    # plots
    'plotly==4.6.0',
    'colorlover',
    # geneditidapp
    'pyramid',
    'pyramid_chameleon',
    'pyramid_jinja2',
    'pyramid_debugtoolbar',
    'pyramid_tm',
    'deform',
    'transaction',
    'zope.sqlalchemy',
    'waitress',
    # amplicount and variantid analysis
    'pyfaidx',
    'biopython',
    'varcode'
]

setup(
    name='geneditid',
    version=version,
    author='GeneditID team',
    description='GeneditID database, WebApp and tools to identify variants and create plate-base visualisation to identify clones with desirable gene edits',
    url='https://geneditid.github.io/',
    license='MIT License',
    packages=find_packages(),
    python_requires='>=3.0, <4',
    include_package_data=True,
    install_requires=install_requires,
    entry_points={
        'paste.app_factory': [
            'main = geneditidapp:main',
        ],
        'console_scripts': [
            'geneditid_create_amplicount_config = scripts.create_ampli_count_conf:main',
            'geneditid_run_amplicount = scripts.run_ampli_count:main',
            'geneditid_run_variantid = scripts.run_variant_id:main',
            'geneditid_add_sample_location = scripts.add_sample_loc:main',
            'geneditid_plot_scores = scripts.plot_scores:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: Python :: 3',
        'Framework :: Pyramid',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],

)

print("""
--------------------------------
GeneditID installation complete!
--------------------------------
For help in running GeneditID, please see the documentation available
at https://geneditid.github.io/
""")
