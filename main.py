from utils import *

from dotenv import load_dotenv; load_dotenv()
openai.api_key = os.getenv("OPENAI_API_KEY")

def main():

    messages = [
        {
        "role": "system",
        "content": "You are an expert computational materials scientist. \nYou are trying to find the EXACT original conditions that were used to make a certain carbon material. You can try varying the number of timesteps and melting temperatures (and cooling rate).\nPredict appropriate code for a melt-quench simulation for this carbon system. \nTry NOT to use an extremely large number of timesteps.\nThe format of the code should be a LAMMPS input file like:\n# Phase PHASE:\nfix integrate all nvt temp T_INITIAL T_FINAL 0.1\nrun TIMESTEPS\nunfix integrate\n"
        },
        {
        "role": "user",
        "content": "EXPERIMENT OUTCOME:\nSOAP Similarity to graphite is 0.1768 \nSOAP Similarity to diamond is 0.9408 \nSOAP Similarity to TARGET is 0.9901\n\nSCRIPT:"
        },
        {
        "role": "assistant",
        "content": "# Phase 1: \nfix integrate all nvt temp 5000 5000 0.1\nrun 50000\nunfix integrate\n\n# Phase 2: \nfix integrate all nvt temp 5000 3000 0.1\nrun 50000 \nunfix integrate\n\n# Phase 3: \nfix integrate all nvt temp 3000 1000 0.1\nrun 50000 \nunfix integrate\n\n# Phase 4:\nfix integrate all nvt temp 1000 300 0.1\nrun 50000 \nunfix integrate\n\n# Phase 5:\nfix integrate all nvt temp 300 300 0.1\nrun 50000\nunfix integrate"
        },
    ]
    QUERY = {
        "role": "user",
        "content": "EXPERIMENT OUTCOME:\nSOAP similarity to graphite is 0.0137\nSOAP similarity to diamond is 0.9820\nSOAP Similarity to TARGET is 1.000\n\nSCRIPT:"
    }

    for i in range(20):
        # format beginning of experiment script
        rundir=abspath(f"run/run_id{i}")
        if not os.path.isdir(rundir): os.makedirs(rundir)
        else:
            # clear existing files and directiories in run directory
            os.system(f'rm -rf {rundir}/*')

        lammps_experiment = format_experiment(rundir)

        # query llm for remaining phases
        RESPONSE, code, oai_response = query_llm(QUERY, messages)
        lammps_experiment += code

        # write experiment to file
        infile = os.path.join(rundir, 'infile')
        with open(infile, 'w') as f:
            f.write(
                lammps_experiment
            )

        # submit and await experiment
        print('Executing experiment...')
        print(RESPONSE)
        print(code)
        sf = '/u/vld/univ5120/VLD/Thesis/llms-for-lammps/lammps/submit.sh'
        # sf = 'lammps/submit.sh'
        schedule(submitfile=sf, infile=infile, rundir=rundir, interval=240)

        # analyse results
        at = read([f for f in glob.glob(os.path.join(rundir, 'dump', 'dump_file.*'))][-1])
        analyses, analysis_nlp = analyse(at)
        # update user response
        messages.append({
            'role': 'user',
            'content': analysis_nlp 
        })
        messages.append(RESPONSE)

        # Logging of query response pairs.
        print('Updated information:\n\n\n', messages)
        with open(os.path.join(rundir,'messages.json'), 'w') as f:
            json.dump(messages, f, indent=4)


if __name__ == '__main__':
    main()


