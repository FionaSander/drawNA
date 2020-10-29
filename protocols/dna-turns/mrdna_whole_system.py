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
#from drawNA.oxdna import Nucleotide

import subprocess
sys.path.insert(0, '/home/fiona/fonso_runs/')
from post_process_copy import *
from functions import initialise_fig
#from transformation import *

FENE_LENGTH = 0.76

#oxDNA_PDB_path = "/Users/fionasander/Desktop/tacoxDNA-python3/src/oxDNA_PDB.py"

def main(length=16, n_strands=10, stapled=5):
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

    # using the two previously created strands create a new strand that
    # we will add to the system
    nucleotides = []
    for strand in strands:
        nucleotides += strand.nucleotides
    strand = Strand(nucleotides=nucleotides)

    # create system and add the final completed strand
    system = System(np.array([100., 100., 100.]))
    system.add_strand(strand)

    actual_doubles = []
    for strand in doubles:
        nucleotides = strand.nucleotides[:stapled]
        actual_doubles.append(Strand(nucleotides=nucleotides))

    system.add_strands(actual_doubles)
    system.write_oxDNA('turns')

if __name__ == '__main__':
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
    parameters['coarse-steps'] = 1e4
    parameters['fine-steps'] = 1e2
    parameters['output-period'] = 1e2

    #parameters['coarse-steps2'] = 1e4
    #parameters['fine-steps2'] = 1e4
    #parameters['output-period2'] = 1e2

    #path = "/home/fiona/fonso_runs" #directory where mrdna simulation output is stored
    #indentify number of simulations carried out
    #simulations = get_simulations(path, 'sim')
    ## Initialize figure with subplots
    #fig = initialise_fig(simulations)
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

    #indentify number of simulations carried out
    #simulations = get_simulations(path, 'sim')

    # Initialize figure with subplots
    fig = initialise_fig(file_list)

    simulation_index = 0 

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
                strand_length = total_length[i]/n_strands
                #file_path = "sim{}_l{}ds{}p{}".format(simulation_index, total_length[i], int(n_strands[j]), stapled[k])# % sim
                #file_path = "sim{}_l{}ds{}p{}".format(simulation_index, total_length[i], stapled[k], int((n_strands*stapled[k])/total_length[i]) )# % sim
                #file_path = "sim{}_l{}ds{}p{}".format(simulation_index, total_length[i], stapled[k], (n_strands*stapled[k])/total_length[i] )# % sim
                file_path = "sim{}_l{}ds{}p{}".format(simulation_index, total_length[i], stapled[k], int(final_percentage*100) )# % sim
                sim_path = file_path + "/sim"
                parameters['fpath'] = file_path
                try:
                    os.mkdir(file_path)
                except FileExistsError:
                    print('{} is being deleted'.format(file_path))
                    shutil.rmtree(file_path)
                    os.mkdir(file_path)

                #Generate oxDNA output files .conf and .top
                #system = main(total_length[i], int(n_strands[j]), stapled[k])
                #system = main(total_length[i], int(n_strands), stapled[k])
                system = main(int(strand_length), int(n_strands), stapled[k])

                shutil.move('./oxdna.turns.top', file_path + '/' + './oxdna.turns.top')
                shutil.move('./oxdna.turns.conf', file_path + '/' + './oxdna.turns.conf')

                #Convert oxDNA output files to .pdb input file format
                subprocess.call(
                    "bash -c \"source ~/.bashrc && module load CUDA && cd {} && python3 /home/fiona/tacoxDNA-python3/src/oxDNA_PDB.py oxdna.turns.top oxdna.turns.conf 35\"".format(file_path), shell=True)

                #mother_system = System(np.array([100., 100., 100.]))

                #i = 0
                #for i in range(4):
                #    i += 1

                #Running mrdna simulation with the .pdb input file
                subprocess.call(
                    "bash -c \"source ~/.bashrc && module load CUDA && cd {} &&".format(parameters['fpath']) + " mrdna --coarse-steps {} --fine-steps {} --output-period {} -d sim oxdna.turns.conf.pdb\"".format(parameters['coarse-steps'],
                                                                                                                                                                                        parameters['fine-steps'],
                                                                                                                                                                                        parameters['output-period']),
                    shell=True)

                    #Convert .pdb intermediate output to oxDNA format (.conf and .top)
            #        subprocess.call(
            #            "bash -c \"source ~/.bashrc && module load CUDA && cd {} && python3 /home/fiona/tacoxDNA-python3/src/PDB_oxDNA.py oxdna.turns.conf-3.pdb 35\"".format(sim_path),
            #            shell=True)   

            #        translation_vector = np.array([10.0, 10.0, 10.0])

            #        strand_system = import_oxDNA(sim_path + "/oxdna.turns.conf-3.pdb.oxdna", sim_path + "/oxdna.turns.conf-3.pdb.top")
            #        #print(strand_system)
            #        new_system = transform(strand_system, translation_vector)
            #        #print(new_system)
            #        new_system.write_oxDNA('turns')

            #        #Convert oxDNA output files to .pdb input file format
            #        subprocess.call(
            #        "bash -c \"source ~/.bashrc && module load CUDA && cd {} && python3 /home/fiona/tacoxDNA-python3/src/oxDNA_PDB.py oxdna.turns.top oxdna.turns.conf 35\"".format(file_path),
            #            shell=True)


            #    mother_system.write_oxDNA('turns')


                #Running mrdna simulation with the .pdb input file from the MOTHER_SYSTEM
               # subprocess.call(
               #     "bash -c \"source ~/.bashrc && module load CUDA && cd {} &&".format(parameters['fpath']) + " mrdna --coarse-steps {} --fine-steps {} --output-period {} -d sim oxdna.turns.conf.pdb\"".format(parameters['coarse-steps2'],
               #                                                                                                                                                                           parameters['fine-steps2'],
               #                                                                                                                                                                           parameters['output-period2']),
               #     shell=True)

                #Running vmd with simulation output to extract .pdb coordinates particles (beads)
                subprocess.call(
                    "bash -c \"source ~/.bashrc && cd {} && vmd -dispdev text -e /home/fiona/fonso_runs/save_all_frames.tcl -args sim\"".format(file_path),shell=True)
                

        

                post_process(file_path, path, file_list, simulation_index, fig)
    
    fig.write_image("/home/fiona/fonso_runs/Rose_diagrams.svg") 
   
