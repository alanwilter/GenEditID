import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

requires = [
    'SQLAlchemy',
]

tests_require = [
    'pytest',
    'pytest-cov',
]

setup(
    name='dnascissors',
    version='1.0',
    description='Gene Editing CRISPR Model',
    classifiers=[
        'Programming Language :: Python',
    ],
    author='',
    author_email='',
    url='',
    keywords='sqlalchemy',
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    extras_require={
        'testing': tests_require,
    },
    install_requires=requires,
)
