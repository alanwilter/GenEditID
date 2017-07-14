from distutils.core import setup

setup(
    name='dnascissors',
    version='1.0',
    author='dnascissors team',
    description='Genome Editing Model',
    long_description=open('README.md').read(),
    license='MIT License',
    packages=['dnascissors', ],
    install_requires=[
        "pyyaml",
        "psycopg2==2.7",
        "sqlalchemy==1.1.6",
        "pandas"
    ],
)
