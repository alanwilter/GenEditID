import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

requires = [
    'pyramid',
    'pyramid_chameleon',
    'pyramid_jinja2',
    'pyramid_debugtoolbar',
    'pyramid_tm',
    'deform',
    'transaction',
    'zope.sqlalchemy',
    'waitress',
    # configuration
    "pyyaml",
    # database
    "psycopg2==2.7",
    'sqlalchemy==1.1.6',
    # for plotting
    'plotly',
    'colorlover',
    # data
    'pandas',
    'xlrd'
]

tests_require = [
    'WebTest >= 1.3.1',  # py3 compat
    'pytest',
    'pytest-cov',
]

setup(
    name='dnascissors-webapp',
    version='1.0',
    author='dnascissors team',
    license='MIT License',
    description='Genome Editing WebApp',
    long_description=open(os.path.join(here, 'README.md')).read(),
    classifiers=[
        'Programming Language :: Python',
        'Framework :: Pyramid',
        'Topic :: Internet :: WWW/HTTP',
        'Topic :: Internet :: WWW/HTTP :: WSGI :: Application',
    ],
        keywords='web pyramid pylons sqlalchemy',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    extras_require={
        'testing': tests_require,
    },
    install_requires=requires,
    entry_points={
        'paste.app_factory': [
            'main = webapp:main',
        ],
        'console_scripts': [
            'initialize_webapp_db = webapp.scripts.initializedb:main',
        ],
    },
)
