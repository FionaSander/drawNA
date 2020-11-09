#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday October 02 2020

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
from drawNA.oxdna import Nucleotide
import subprocess
sys.path.insert(0, '/home/fiona/fonso_runs/')
from post_process_copy import *
from functions import initialise_fig
from transformation import *
from copy import copy

FENE_LENGTH = 0.76

#Generation of initial strand configuration
def generate_system(length=16, n_strands=10, stapled=5):
    # generate a strand
    strands = []
    doubles = []
    strand, double = generate_helix(n=length, double=True)
    strands.append(strand.copy())
    doubles.append(double.copy())

    for i in range(n_strands-1):

        last_nuc = strands[-1].nucleotides[-1]
        direction = -last_nuc._a3
        a1 = -last_nuc._a1

        # ensure the backbone position is FENE_LENGTH away from
        # the backbone position of the previous nucleotide
        start = last_nuc.pos_back + (FENE_LENGTH - POS_BACK) * a1

        # generate strand above that's going in opposite direction
        strand, double = generate_helix(
            n=length,
            start_position=start,
            direction=direction,
            a1=a1,
            double=True,
        )
        strands.append(strand)
        doubles.append(double)

    # using the two previously created strands create a new strand that is added to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides
    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    main_system = System(np.array([80.0, 80.0, 80.0]))
    main_system.add_strand(strand)

    actual_doubles = []
    for strand in doubles:
        nucleotides = strand.nucleotides[:stapled]
        actual_doubles.append(Strand(nucleotides=nucleotides))

    main_system.add_strands(actual_doubles)
    
    #main_system.write_oxDNA('turns')

    return main_system


# conversion of drawNA output to mrdna input parameters
#if __name__ == '__main__':
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-l', '--length', type=int, nargs=3, required=True)
parser.add_argument('-ds', '--double-stranded', type=int, nargs=3, required=True)
parser.add_argument('-p', '--double-stranded-percentage', type=int, nargs=3, required=True)
args = parser.parse_args()

total_length = np.arange(args.length[0], args.length[2]+1, args.length[1])
percentage = np.arange(args.double_stranded_percentage[0], args.double_stranded_percentage[2]+1, args.double_stranded_percentage[1])
stapled = np.arange(args.double_stranded[0], args.double_stranded[2]+1, args.double_stranded[1])

simulation_index = 0  
parameters = dict()
parameters['coarse-steps'] = 1e5
parameters['fine-steps'] = 1e5
parameters['output-period'] = 1e2

parameters['coarse-steps2'] = 1e5
parameters['fine-steps2'] = 1e5
parameters['output-period2'] = 1e2


file_list = []
for i in range (len(total_length)):
    #strand_length = total_length[i]*(1.0/percentage)
    #n_strands = total_length[i]/strand_length
    for j in range (len(percentage)):
        for k in range (len(stapled)):
            bases_needed = total_length[i] * (1/(100.0/percentage[j]))
            n_strands = np.round(bases_needed/stapled[k])

            #Create simulation output folder
            simulation_index = simulation_index + 1
            final_percentage = round(n_strands*stapled[k])/total_length[i]
            #print(final_percentage)
            #file_list.append("sim{}_l{}ds{}p{}".format(simulation_index, total_length[i], int(n_strands[j]), stapled[k]))
            #file_list.append("sim{}_l{}ds{}p{}".format(simulation_index, total_length[i], stapled[k], int((n_strands*stapled[k])/total_length[i]) ))# % sim
            file_list.append("sim{}_l{}ds{}p{}".format(simulation_index, total_length[i], stapled[k], int(final_percentage*100)) )# % sim
# define the different directories here
path = "/home/fiona/fonso_runs" #directory where mrdna simulation output is stored

# Initialize figure with subplots
fig = initialise_fig(file_list)

simulation_index = 0 

for l in range (len(total_length)):
    for j in range (len(percentage)):
        for k in range (len(stapled)):
            bases_needed = total_length[l] * (1/(100.0/percentage[j]))
            n_strands = np.round(bases_needed/stapled[k])

            #Create simulation output folder
            simulation_index = simulation_index + 1
            final_percentage = round(n_strands*stapled[k])/total_length[l]
            strand_length = total_length[l]/n_strands
            file_path = "sim{}_l{}ds{}p{}".format(simulation_index, total_length[l], stapled[k], int(final_percentage*100) )# % sim
            sim_path = file_path + "/sim"
            parameters['fpath'] = file_path
            try:
                os.mkdir(file_path)
            except FileExistsError:
                print('{} is being deleted'.format(file_path))
                shutil.rmtree(file_path)
                os.mkdir(file_path)

            #Generate oxDNA output files .conf and .top
            main_system = generate_system(int(strand_length), int(n_strands), stapled[k])

    
            # Definitions
            DELTA = np.array([0., 0., 10.])
            MAIN = 'main'
            TEMPLATE = MAIN + '.{}'
            N = 3
            SIMULATIONS = range(1, N)


            # create clone_system
            clone_system = copy(main_system)
            clone_name = TEMPLATE.format(0)
            clone_pdb = f'oxdna.{clone_name}.conf.pdb'

            # write clones_system to oxDNA
            clone_system.write_oxDNA(clone_name)

            #move .conf and .top files to correct folder
            shutil.move('./oxdna.{}.top'.format(clone_name), file_path + '/' + './oxdna.{}.top'.format(clone_name))
            shutil.move('./oxdna.{}.conf'.format(clone_name), file_path + '/' + './oxdna.{}.conf'.format(clone_name))

            #Convert oxDNA output files to .pdb input file format
            subprocess.call(
            "bash -c \"source ~/.bashrc && module load CUDA && cd {} && python3 /home/fiona/tacoxDNA-python3/src/oxDNA_PDB.py oxdna.{}.top oxdna.{}.conf 35\"".format(file_path, clone_name, clone_name),
                shell=True)

            os.remove(file_path + '/oxdna.{}.conf'.format(clone_name))
            os.remove(file_path + '/oxdna.{}.top'.format(clone_name))

            for i in SIMULATIONS:

                new_name = TEMPLATE.format(i)
                new_pdb = f'oxdna.{new_name}.conf.pdb'

                # Running mrdna simulation
                run_simulation(clone_pdb, new_pdb, file_path, parameters) 
                shutil.rmtree(sim_path)

                # convert pdb to oxDNA and import into drawNA
                clone_system = read_pdb(clone_pdb, sim_path, file_path)

                # remove not needed oxDNA and pdb files
                os.remove(file_path + '/{}.top'.format(clone_pdb))
                os.remove(file_path + '/{}.conf'.format(clone_pdb))
                os.remove(file_path + '/{}'.format(clone_pdb))

                # overwrite clone_name and clone_pdb with new system
                clone_name = new_name
                clone_pdb = new_pdb

                # translate the new system
                translate_system(clone_system, DELTA)
                DELTA[2] = DELTA[2] + 10.0

                # clone the strands and add them to the main system
                clone_strands = clone_system.strands
                main_system.add_strands(clone_strands)
                
                # write oxDNA system
                clone_system.write_oxDNA(clone_name)

                # move oxDNA files to correct folder
                shutil.move('./oxdna.{}.top'.format(new_name), file_path + '/' + './oxdna.{}.top'.format(new_name))
                shutil.move('./oxdna.{}.conf'.format(new_name), file_path + '/' + './oxdna.{}.conf'.format(new_name))

                # convert oxDNA to pdb
                convert_to_pdb(clone_name, clone_pdb, file_path, i)

                # remove not needed oxDNA files
                os.remove(file_path + '/oxdna.{}.top'.format(new_name))
                os.remove(file_path + '/oxdna.{}.conf'.format(new_name))

            # remove not needed pdb files
            os.remove(file_path + '/oxdna.{}.conf.pdb'.format(new_name))

            # write main system to oxDNA
            main_system.write_oxDNA('MAIN')

            # move oxDNA files to correct folder
            shutil.move('./oxdna.MAIN.top', file_path + '/oxdna.MAIN.top')
            shutil.move('./oxdna.MAIN.conf', file_path + '/oxdna.MAIN.conf')

            #Convert oxDNA output files to .pdb input file format
            subprocess.call(
            "bash -c \"source ~/.bashrc && module load CUDA && cd {} && python3 /home/fiona/tacoxDNA-python3/src/oxDNA_PDB.py oxdna.MAIN.top oxdna.MAIN.conf 35\"".format(file_path),
                shell=True)

            #Running mrdna simulation with the .pdb input file from the MOTHER_SYSTEM
            subprocess.call(
                "bash -c \"source ~/.bashrc && module load CUDA && cd {} &&".format(parameters['fpath']) + " mrdna --coarse-steps {} --fine-steps {} --output-period {} -d sim oxdna.MAIN.conf.pdb\"".format(parameters['coarse-steps2'],
                                                                                                                                                                                        parameters['fine-steps2'],
                                                                                                                                                                                        parameters['output-period2']),
                shell=True)

            #Running vmd with simulation output to extract .pdb coordinates particles (beads)
            subprocess.call(
                "bash -c \"source ~/.bashrc && cd {} && vmd -dispdev text -e /home/fiona/fonso_runs/save_all_frames.tcl -args sim\"".format(file_path),shell=True)
            
            post_process(file_path, path, file_list, simulation_index, fig)

fig.write_image("/home/fiona/fonso_runs/Rose_diagrams.svg") 

