
# Base
import os, sys, glob
from os.path import abspath

# Chemistry
import ase
from ase.io import read, write

# ML / Linalg
import numpy as np


def format_query(past_queries):

    sys_msg = """You are an expert computational materials scientist. 
You are trying to develop an experiment to crystallise carbon into diamond."""

    q = """Predict appropriate code for a melt-quench simulation of this system. 
The format of the code should be a LAMMPS input file.
A 

# EXAMPLE 1
RESULT
SOAP similarity to diamond is 0.9
SOAP similarity to graphite is 0.1

SCRIPT:
# Phase 1: Melt
fix integrate all nvt temp 9000 9000 0.1
run 100000
unfix integrate

# Phase 2: Quench
fix integrate all nvt temp 9000 3000 0.1
run 100000 
unfix integrate

# Phase 3: Warm
fix integrate all nvt temp 3000 500 0.1
run 100000 
unfix integrate

# Phase 4: Anneal
fix integrate all nvt temp 500 500 0.1
run 100000 
unfix integrate

# Phase 5: Cool
fix integrate all nvt temp 500 300 0.1
run 100000
unfix integrate

# EXAMPLE 2 
RESULT:
SOAP similarity to diamond is 1.0
SOAP similarity to graphite is 0.0

SCRIPT:
"""
    
    ##########
    # Include a number of past experiments
    q += ""  #  TODO: format past queries
    ##########
    

    return q
    
def infer_experiment_phase(query):
    # Use LLM to infer the next phase of the MD
    ...

def format_experiment(query):

    # MQ params
    # timesteps in ps
    melt=5
    quench=5
    warm=50
    anneal=100
    cool=50

    # temperatures in K
    melt_T=9000
    quench_T=3000

    structure_file = abspath("structures/example.extxyz") 
    dump_file = abspath("dump_file")
    pot = abspath("data/carbon.xml")
    rand_seed = 42  

    # Step 1: Declare necessary variables
    quench = 100
    warm = 100
    anneal = 100 
    cool = 100
    melt = 100

    melt_T = 9000
    quench_T = 3000
    anneal_T = 500
    cool_T = 300


    # Step 2: Create LAMMPS input file string
    lammps_input = f"""
# variable melt_steps equal {melt*1000}
# variable quench_steps equal {quench*1000}
# variable warm_steps equal {warm*1000}
# variable anneal_steps equal {anneal*1000}
# variable cool_steps equal {cool*1000}

# variable Nfreq equal 100
# variable Nevery equal 100
# variable Nrepeat equal {100/(warm*1000)}
# variable Ndump equal 1000

print "Melt  : {melt*1000}fs"
print "Quench: {quench*1000}fs"
print "Warm  : {warm*1000}fs"
print "Anneal: {anneal*1000}fs"
print "Cool  : {cool*1000}fs"

log ./log.dat append

units metal     # mass = g/mole, distance = Å, time = ps, energy = eV, force = eV/Å
atom_style atomic
read_data {structure_file}
mass * 12.011   # all atoms are Carbon and therefore have mass ~12

pair_style quip
pair_coeff * * {pot} "" 6

neigh_modify every 1 delay 0 check yes
timestep 0.001

variable nAtoms equal atoms

fix removeMomentum all momentum 1 linear 1 1 1

compute T all temp
fix TempAve all ave/time {warm*1000} {100/(warm*1000)} 100 c_T

variable P equal press
fix PressAve all ave/time {warm*1000} {100/(warm*1000)} 100 v_P

variable v equal vol
fix vAve all ave/time {warm*1000} {100/(warm*1000)} 100 v_v

compute PE all pe pair

variable PE_Atom equal c_PE/v_nAtoms

fix PEAve_Atom all ave/time {warm*1000} {100/(warm*1000)} 100 v_PE_Atom

compute MSD all msd

thermo_style custom step cpu temp f_TempAve press f_PressAve f_PEAve_Atom vol f_vAve c_MSD[4]
thermo_modify flush yes
thermo {warm*1000}

dump traj all atom {warm*1000} {dump_file}.dump.*.dat

velocity all create {melt_T} {rand_seed}

run 0

# Trajectory Phases
# Phase 1: Melt
fix integrate all nvt temp {melt_T} {melt_T} 0.1
run {melt*1000}
unfix integrate

# Phase 2: Quench
fix integrate all nvt temp {melt_T} {quench_T} 0.1
run {quench*1000} 
unfix integrate

# Phase 3: Warm
fix integrate all nvt temp {quench_T} {anneal_T} 0.1
run {warm*1000} 
unfix integrate

# Phase 4: Anneal
fix integrate all nvt temp {anneal_T} {anneal_T} 0.1
run {anneal*1000} 
unfix integrate

# Phase 5: Cool
fix integrate all nvt temp {anneal_T} {cool_T} 0.1
run {cool*1000}
unfix integrate

"""
    return lammps_input


def main():
    past_queries = [
        {
            'lammps_phase': '',
            'metrics': {
                'soap_similarity': None,
            }
        }
    ]
    while True:
        # (1) Get updated query
        query = format_query(query_result)
        # (2) Infer experiment from query
        experiment_string = infer_experiment(query)
        # (3) Run experiment
        print('submitting experiment: \n\n\n', experiment_string)
        with open('lammps/experiment.in', 'w') as f:
            f.write(experiment_string)
        os.system('qsub submit.sh lammps/experiment.in')
        # Check SGE queue for job completion
        atoms = read('')

        # (4) Compute metrics
        # (5) Update db
        break



if __name__ == '__main__':
    # few tests
    # x = infer_experiment(query)
    s = format_experiment(None)
    print(s)
    infile_path = '/home/jbutch/Projects/Thesis/llms-for-lammps/lammps/experiment.in'
    with open(infile_path, 'w') as f:
        f.write(s)

    # main()

