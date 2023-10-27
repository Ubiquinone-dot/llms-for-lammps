
# Base
import os, sys, glob
from os.path import abspath

# Chemistry
import ase
from ase.io import read, write

# ML / Linalg
import numpy as np

# OAI
import openai
from dotenv import load_dotenv; load_dotenv()
openai.api_key = os.getenv("OPENAI_API_KEY")

def query_llm(messages=[]):

    response = openai.ChatCompletion.create(
    model="gpt-4",
    messages=messages,
        temperature=0.1,
        max_tokens=256,
        top_p=1,
        frequency_penalty=0,
        presence_penalty=0
    )

    message = response.choices[0].message
    messages.append(message)
    code = message['content']
    print(response)
    print(messages)
    print(code)
    return code, messages, response
    

def format_experiment():

    structure_file = abspath("structures/example.extxyz") 
    dump_file = abspath("dump_file")  # no extension
    log_file = abspath("log.dat")
    pot = abspath("data/carbon.xml")
    rand_seed = 42  

    ncomp = 100000  # (1000*nps) in fs
    compfreq = 100/(ncomp*1000)
    T_0 = 3000  # initial temperature for velcoity creation

    lammps_input = f"""
log {log_file} append

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
variable P equal press
variable v equal vol
variable PE_Atom equal c_PE/v_nAtoms

fix TempAve all ave/time {ncomp} {compfreq} 100 c_T
fix PressAve all ave/time {ncomp} {compfreq} 100 v_P
fix vAve all ave/time {ncomp} {compfreq} 100 v_v
fix PEAve_Atom all ave/time {ncomp} {compfreq} 100 v_PE_Atom

compute PE all pe pair
compute MSD all msd

thermo_style custom step cpu temp f_TempAve press f_PressAve f_PEAve_Atom vol f_vAve c_MSD[4]
thermo_modify flush yes
thermo {ncomp}

dump traj all atom {ncomp} {dump_file}.dump.*.dat

velocity all create {T_0} {rand_seed}

run 0

# PHASES:
"""
    # When a full experiment is requested, add some phases to the experiment
#     if FULL: lammps_input += f"""# Prebuilt Phases
# # Phase 1: Melt
# fix integrate all nvt temp {melt_T} {melt_T} 0.1
# run {melt*1000}
# unfix integrate

# # Phase 2: Quench
# fix integrate all nvt temp {melt_T} {quench_T} 0.1
# run {quench*1000} 
# unfix integrate

# # Phase 3: Warm
# fix integrate all nvt temp {quench_T} {anneal_T} 0.1
# run {warm*1000} 
# unfix integrate

# # Phase 4: Anneal
# fix integrate all nvt temp {anneal_T} {anneal_T} 0.1
# run {anneal*1000} 
# unfix integrate

# # Phase 5: Cool
# fix integrate all nvt temp {anneal_T} {cool_T} 0.1
# run {cool*1000}
# unfix integrate

# """
    return lammps_input


def main():

    messages = [
        {
        "role": "system",
        "content": "You are an expert computational materials scientist. \nYou are trying to develop an experiment to crystallise carbon into diamond.\nPredict appropriate code for a melt-quench simulation of this system. \nTry NOT to use an extremely large number of timesteps.\nThe format of the code should be a LAMMPS input file like:\n# Phase PHASE:\nfix integrate all nvt temp T_INITIAL T_FINAL 0.1\nrun TIMESTEPS\nunfix integrate\n"
        },
        {
        "role": "user",
        "content": "RESULT:\nSOAP similarity to diamond is 0.1\nSOAP similarity to graphite is 0.9\n\nSCRIPT:"
        },
        {
        "role": "assistant",
        "content": "# Phase 1: \nfix integrate all nvt temp 9000 9000 0.1\nrun 100000\nunfix integrate\n\n# Phase 2: \nfix integrate all nvt temp 9000 3000 0.1\nrun 100000 \nunfix integrate\n\n# Phase 3: \nfix integrate all nvt temp 3000 500 0.1\nrun 100000 \nunfix integrate\n\n# Phase 4:\nfix integrate all nvt temp 500 500 0.1\nrun 100000 \nunfix integrate\n\n# Phase 5:\nfix integrate all nvt temp 500 300 0.1\nrun 100000\nunfix integrate"
        },
        {
        "role": "user",
        "content": "RESULT:\nSOAP similarity to diamond is 1.0\nSOAP similarity to graphite is 0.0\n\nSCRIPT:"
        },
    ]
    while True:
        # query llm
        code, messages, response = query_llm(messages)

        # submit experiment

        # wait

        # analyse results and update system response

        exit()

if __name__ == '__main__':
    # # TEST FORMATTING
    # s = format_experiment()
    # print(s)
    # infile_path = '/home/jbutch/Projects/Thesis/llms-for-lammps/lammps/experiment_human.in'
    # with open(infile_path, 'w') as f:
    #     f.write(s)

    phases_str = query_llm()
    print(phases_str)

    s = format_experiment(FULL=False)
    s += phases_str  # append phases to experiment
    infile_path = '/home/jbutch/Projects/Thesis/llms-for-lammps/lammps/experiment_ai.in'
    with open(infile_path, 'w') as f:
        f.write(s)
    print(s)

    code, messages, response
