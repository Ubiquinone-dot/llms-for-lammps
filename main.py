from utils import *


def main():

    messages = [
        {
        "role": "system",
        "content": "You are an expert computational materials scientist. \nYou are trying to find the EXACT original conditions that were used to make a material. You can try varying the number of timesteps and melting temperatures (and thus cooling rate).\nPredict appropriate code for a melt-quench simulation for this carbon system. \nTry NOT to use an extremely large number of timesteps.\nThe format of the code should be a LAMMPS input file like:\n# Phase PHASE:\nfix integrate all nvt temp T_INITIAL T_FINAL 0.1\nrun TIMESTEPS\nunfix integrate\n"
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
        "content": "RESULT:\nSOAP similarity to diamond is 1.0\nSOAP similarity to graphite is 0.03\n\nSCRIPT:"
        },
    ]

    for i in range(2):
        # format beginning of experiment script
        rundir=abspath(f"run/run_id{run_id}")

        if not os.path.isdir(rundir): os.makedirs(rundir)
        lammps_experiment = format_experiment(rundir)

        # query llm for remaining phases
        code, messages, response = query_llm(messages)
        lammps_experiment += code

        # write experiment to file
        infile = os.path.join(rundir, 'infile')
        with open(infile, 'w') as f:
            f.write(
                lammps_experiment
            )

        # submit and await experiment
        print('Executing experiment...')
        print(response)
        print(code)
        schedule(submitfile='lammps/submit.sh', infile=infile, rundir=rundir, interval=240)

        # analyse results
        at = read([f for f in glob.glob(os.path.join(rundir, 'dump', 'dump_file.*'))][-1])
        analyses, analysis_nlp = analyse(at)
        print(sims_str)

        # update user response
        messages.append({
            'role': 'user',
            'content': analysis_nlp 
        })
        print(messages)


if __name__ == '__main__':
    main()


# # TEST FORMATTING
# s = format_experiment()
# print(s)
# infile_path = '/home/jbutch/Projects/Thesis/llms-for-lammps/lammps/experiment_human.in'
# with open(infile_path, 'w') as f:
#     f.write(s)

# phases_str = query_llm()
# print(phases_str)

# s = format_experiment(FULL=False)
# s += phases_str  # append phases to experiment
# infile_path = '/home/jbutch/Projects/Thesis/llms-for-lammps/lammps/experiment_ai.in'
# with open(infile_path, 'w') as f:
#     f.write(s)
# print(s)

# code, messages, response
