#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Monday October 26 2020

@author: fionasander
"""
import os
import sys
print(os.path.dirname(__file__) + '/../../')
sys.path.insert(0, os.path.dirname(__file__) + '/../../')
import numpy as np
import shutil, os
from drawNA.oxdna.strand import generate_helix, POS_BACK, Strand
from drawNA.oxdna import System
from drawNA.readers import OXDNAReader
from drawNA.oxdna import Nucleotide
import subprocess
sys.path.insert(0, '/home/fiona/fonso_runs/')
from post_process_copy import *
from functions import initialise_fig
from copy import deepcopy


def import_oxDNA(conf, top):
    """
    import .conf and .top files to drawNA system
    """
    reader = OXDNAReader([conf, top])
    strand_system = reader.system #or reader.strand??

    return strand_system


def add_strand_system_to_mother_system(mother_system, strand_system):
    """
    add an imported strand from mrdna to the "mother" system
    """
    mother_system.add_strand(strand_system)

    return mother_system


#def transform(strand_system: drawNA.oxdna.Nucleotide, translation_vector: np.ndarray):
def transform(strand_system, translation_vector: np.ndarray):
    """
    shift all coordinates of a system by a translation vector
    """
    new_system = deepcopy(strand_system)

    for nt in new_system._nucleotides:
        nt.pos_com += translation_vector

    return new_system




