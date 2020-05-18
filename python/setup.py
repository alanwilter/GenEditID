from setuptools import setup, find_packages

version = '0.2'

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
    # test
    'pytest',
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
            'geneditid_run_ampli_find = geneditidtools.run_ampli_find:main',
            'geneditid_run_amplicount = geneditidtools.run_ampli_count:main',
            'geneditid_run_variantid = geneditidtools.run_ampli_plots:main',
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
