# Util functions for main.py
# Base
import os, sys, glob, json
from os.path import abspath
import subprocess
import time, datetime
import re

# Chemistry
import ase
from ase.io import read, write
from dscribe.descriptors import SOAP
from ase.io import read, write
from ase import build

# ML / Linalg
import numpy as np

# OAI
import openai
from dotenv import load_dotenv; load_dotenv()
openai.api_key = os.getenv("OPENAI_API_KEY")


def query_llm(QUERY, messages=[]):

    response = openai.ChatCompletion.create(
    model="gpt-4",
    messages=messages+[QUERY],
        temperature=0.1,
        max_tokens=256,
        top_p=1,
        frequency_penalty=0,
        presence_penalty=0
    )

    RESPONSE = response.choices[0].message
    code = RESPONSE['content']

    return RESPONSE, code, response
    

def format_experiment(rundir=None):
    if rundir is None:
        rundir = abspath("run/test_run")
    
    structure_file = abspath("structures/example.data") 
    dump_file = abspath(os.path.join(rundir, "dump/dump_file"))  # no extension
    log_file = abspath(os.path.join(rundir, "log.dat"))
    pot_file = abspath("data/carbon.xml")

    if not os.path.exists(rundir): 
        os.makedirs(rundir)
    if not os.path.exists(os.path.join(rundir, "dump")): 
        os.makedirs(os.path.join(rundir, "dump"))
    
    print('formatting experiment for rundir:\n', rundir)

    nevery = 100                # how many timesteps between each fix 
    nfreq = 100                 # how many timesteps between each dump
    nrepeat = int(nfreq/nevery) # how many times to repeat the fix

    rand_seed = 42
    T_0 = 3000  # initial temperature for velcoity creation

    lammps_input = f"""
log {log_file} append

units metal     # mass = g/mole, distance = Å, time = ps, energy = eV, force = eV/Å
atom_style atomic
read_data {structure_file}
mass * 12.011   # all atoms are Carbon and therefore have mass ~12

pair_style quip
pair_coeff * * {pot_file} "" 6

neigh_modify every 1 delay 0 check yes
timestep 0.001

variable nAtoms equal atoms

fix removeMomentum all momentum 1 linear 1 1 1

compute T all temp
variable P equal press
variable v equal vol
variable PE_Atom equal c_PE/v_nAtoms

fix TempAve all ave/time {nevery} {nrepeat} {nfreq} c_T
fix PressAve all ave/time {nevery} {nrepeat} {nfreq} v_P
fix vAve all ave/time {nevery} {nrepeat} {nfreq} v_v
fix PEAve_Atom all ave/time {nevery} {nrepeat} {nfreq} v_PE_Atom

compute PE all pe pair
compute MSD all msd

thermo_style custom step cpu temp f_TempAve press f_PressAve f_PEAve_Atom vol f_vAve c_MSD[4]
thermo_modify flush yes
thermo {nevery}

dump traj all atom {nevery} {dump_file}.dump.*.dat

velocity all create {T_0} {rand_seed}

run 0

# (AI) PHASES: 

"""
    return lammps_input



def schedule(submitfile, infile, rundir, interval=30):
    time.sleep(0.1)
    
    # Submit the job and capture the output
    result = subprocess.run(['qsub', f'{submitfile}', 
        f'{infile}', f'{rundir}'  # arguments to the script $1, $2.
    ], stdout=subprocess.PIPE, text=True, cwd=rundir)
    output = result.stdout
    print('Submission output:', output)

    # Extract the job ID from the output
    job_id_match = re.search(r'(\d+)', output)
    if job_id_match is None:
        print('Failed to submit job or parse job ID')
        exit(1)

    job_id = job_id_match.group(1)
    print('Job ID:', job_id)

    while True:
        # Query the job status
        result = subprocess.run(['qstat', '-j', job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            print(f'{datetime.datetime.now().strftime("%H:%M")}: Job completed or failed')
            break
        
        # TODO: Parse the job status from the output  
        print(f'{datetime.datetime.now().strftime("%H:%M")}: '+f'Job status: {None}')
        # status_match = re.search(r'job_state: (\w+)', result.stdout)
        # if status_match:
        #     status_code = status_match.group(1)
        #     status_description = job_status_codes.get(status_code, 'unknown status')
        #     print(f'{datetime.datetime.now().strftime("%H:%M")}: Job status: {status_description}')
        # else:
        #     print(f'{datetime.datetime.now().strftime("%H:%M")}: Failed to parse job status')

        time.sleep(interval)


def analyse(at, rattle=True):  # written explicity for carbon

    def kernel(soaps1, soaps2, zeta=2):
        kernel_mat = np.dot(soaps1, soaps2)**zeta

        norm1 = np.linalg.norm(soaps1, axis=-1)**zeta
        norm2 = np.linalg.norm(soaps2, axis=-1)**zeta
        prod = norm1 * norm2
        kernel_mat = kernel_mat / prod

        return kernel_mat
        
    # dia_at = read('structures/diamond.xyz')
    # # gra_at = read('structures/graphite.xyz')
    dia_at = build.bulk('C', 'diamond', a=3.567)
    gra_at = build.graphene('C2', a=2.46)
    
    at_original = read('structures/example.extxyz')

    # double check for lammps file loading errors:
    at.set_atomic_numbers([6 for i in range(len(at))])

    if rattle:
        # Supercell and move a little to diversity soaps a little.
        dia_supercell = dia_at * (2, 2, 2)
        gra_supercell = gra_at * (3, 3, 1)

        # Define the magnitude of random perturbation
        perturbation_magnitude = 0.01  # e.g., 0.01 Angstrom

        # Apply random perturbations to atomic positions
        dia_supercell.positions += (np.random.rand(*dia_supercell.positions.shape)) * perturbation_magnitude
        gra_supercell.positions += (np.random.rand(*gra_supercell.positions.shape)) * perturbation_magnitude

    soap_gen = SOAP(species=['C'], periodic=True, r_cut=4.8, n_max=8, l_max=6)
    at_soap = soap_gen.create(at)
    dia_soap = soap_gen.create(dia_at)
    at_original_soap = soap_gen.create(at_original)

    soap_gen = SOAP(species=['C'], r_cut=4.8, n_max=8, l_max=6)
    gra_soap = soap_gen.create(gra_at)
    
    # take average soap vector for each structure
    at_soap = np.mean(at_soap, axis=0)
    dia_soap = np.mean(dia_soap, axis=0)
    gra_soap = np.mean(gra_soap, axis=0)
    at_original_soap = np.mean(at_original_soap, axis=0)

    # krnl = kernel(soaps)
    sims = {
        'graphite': kernel(at_soap, gra_soap),  # could leave out?
        'diamond': kernel(at_soap, dia_soap),
        'TARGET': kernel(at_soap, at_original_soap)
    }
    # make floating point precision to 4 decimal places
    sims = { k: round(v, 4) for k,v in sims.items() }

    sims_str = '\n'.join([ f'SOAP Similarity to {k} is {v} ' for k,v in sims.items() ])  # SOAP similarity to diamond is {similarity_to_diamond}\nSOAP similarity to graphite is {similarity_to_graphite}
    sims_str = 'EXPERIMENT OUTCOME:\n' + sims_str + '\nSCRIPT:'
    
    return sims, sims_str