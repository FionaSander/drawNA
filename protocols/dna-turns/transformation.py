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
from copy import deepcopy, copy


def import_oxDNA(conf, top):
    """
    import .conf and .top files to drawNA system
    """
    reader = OXDNAReader([conf, top])
    strand_system = reader.system #or reader.strands??

    strand_system_strand = strand_system._strands # how do i extract strand from strand_system?

    return strand_system#, strand_system_strand


def add_strand_to_mother_system(mother_system, strand_system_strand):
    """
    add an imported strand from mrdna to the "mother" system
    """
    mother_system.add_strands(strand_system_strand)

    return mother_system


#def transform(strand_system: drawNA.oxdna.Nucleotide, translation_vector: np.ndarray):
def transform(strand_system, translation_vector: np.ndarray):
    """
    shift all coordinates of a system by a translation vector
    """
    new_system = copy(strand_system)

    for nt in new_system.nucleotides:#_nucleotides: not sure if this is the right way to do it but lets see...
        nt.pos_com += translation_vector

    return new_system

def read_pdb(clone_pdb, sim_path, file_path):
    """
    Convert .pdb output to oxDNA format (.conf and .top) and import it into drawNA
    """
    #import explicitly the atomic resolution output file (-3.pdb) - main.0.conf-3.pdb (this is how the files are created in mrdna)
    #clone_pdb = clone_pdb[:-5] + '-3' + clone_pdb[-5:]

    #shutil.move('./{}'.format(clone_pdb), file_path + '/' + './{}'.format(clone_pdb))
    

    #conversion to oxDNA
    subprocess.call(
        "bash -c \"source ~/.bashrc && module load CUDA && cd {} && python3 /home/fiona/tacoxDNA-python3/src/PDB_oxDNA.py {} 35\"".format(file_path, clone_pdb),
        shell=True)   


    #import oxDNA system into drawNA into a strand_system
    #'oxdna.{}.conf.pdb'.format(TEMPLATE.format(i+1))
    clone_system = import_oxDNA( file_path + "/{}.conf".format(clone_pdb), file_path + "/{}.top".format(clone_pdb) )
    #clone_system = import_oxDNA( file_path + "/oxdna.main.{}.conf".format((i)), file_path + "/oxdna.main.{}.top".format((i)) )

    return clone_system

def convert_to_pdb(clone_name, clone_pdb, file_path, i):
    """
    convert drawNA system to PDB file format
    """
    #Convert oxDNA output files to .pdb input file format
    subprocess.call(
    "bash -c \"source ~/.bashrc && module load CUDA && cd {} && python3 /home/fiona/tacoxDNA-python3/src/oxDNA_PDB.py oxdna.{}.top oxdna.{}.conf 35\"".format(file_path, clone_name, clone_name),
        shell=True)
    
    #return clone_pdb for the next loop
    #clone_pdb = "{}{}.conf.pdb".format(file_path, clone_pdb)
    #oxdna.main.1.conf.pdb
    
    #clone_pdb = 'oxdna.main.{}.conf.pdb'.format(TEMPLATE,(i+1))
    clone_pdb = 'oxdna.main.{}.conf.pdb'.format((i))

def translate_system(system, translation_vector: np.ndarray):
    """
    shift all coordinates of a system by a translation vector
    """

    for nt in system.nucleotides:#_nucleotides: not sure if this is the right way to do it but lets see...
        nt.pos_com += translation_vector

    return system


def run_simulation(clone_pdb, new_pdb, file_path, parameters):
    """
    run mrdna simulation
    """
    #Running mrdna simulation with the .pdb input file on one strand
    subprocess.call(
        "bash -c \"source ~/.bashrc && module load CUDA && cd {} &&".format(parameters['fpath']) + " mrdna --coarse-steps {} --fine-steps {} --output-period {} -d sim {}\"".format(parameters['coarse-steps'],
                                                                                                                                                                            parameters['fine-steps'],
                                                                                                                                                                            parameters['output-period'], clone_pdb),
        shell=True)

    
    new_pdb = file_path + clone_pdb[:-5] + '-3' + clone_pdb[-5:] #new_pdb contains the atomisitic pdb from the mrdna simulation output

    #remove /sim folder (not needed anymore for sub systems)


    return new_pdb


