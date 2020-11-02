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
    strand_system = reader.system 

    strand_system_strand = strand_system._strands # is this how i extract strands from strand_system?

    return strand_system, strand_system_strand


def add_strand_to_mother_system(mother_system, strand_system_strand):
    """
    add an imported strand from mrdna to the "mother" system
    """
    mother_system.add_strands(strand_system_strand)

    return mother_system


def transform(strand_system, translation_vector: np.ndarray):
    """
    shift all coordinates of a system by a translation vector
    """
    new_system = deepcopy(strand_system)

    for nt in new_system.nucleotides:#_nucleotides: not sure if this is the right way to do it but lets see...
        nt.pos_com += translation_vector

    return new_system




