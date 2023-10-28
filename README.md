# llms-for-lammps
Uses llms to design experiments to submit to lammps. Runs repeatedly to get toward some structure.


## More details
See the messages variable in main to see what context the LLM is provided with to begin with, alternatively see the playground API [here](https://platform.openai.com/playground/p/xL0CfjMc6pTnUbk5WIjEz4mj?model=gpt-4). The basic idea is that we've got some structure which was made by a melt-quench simulation, and we'll ask the LLM to try and infer the experiment (I've forgotten how it was made so). Obviously it's not given much context but we feed back the results to it and so maybe it can manage to infer the best parameters after a few datapoints.

See: 
- `main.py` for 
- `utils.py` for; (i) formatting the details of the experiment (linking together all the necessary components of LAMMPS) (ii) functionality to submit a job to a SGE compute cluster (Sun-grid-engine) and await the result and (ii) analysis

## Todo's
You can expand the analysis to feed back better information to the LLM

You can change the way we get the LLM to infer the parameters of the simulation to lower the number of tokens required, lowering the cost of running the simulations continuously and avoiding context limits.

You can give additional parameters to the LLM to infer to make things more impressive; change density, size, number of atoms etc in original structure.

