"""
Python3 file of the genome-editing project
Created by Anne Pajon @pajanne on 08/03/2017
"""

import yaml
import os

yml_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'crispr.yml')

with open(yml_filepath, 'r') as yml_file:
    cfg = yaml.load(yml_file, Loader=yaml.FullLoader)
