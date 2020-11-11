"""
This is the master file, which you should use to set up and run the simulations.
You may define functions or classes as necessary

For an input sequence, do the following (see project page for details):
	1. load sequence from fasta file and initialize protein into extended configuration
	2. Run a series of simulations (n=5 or 10):
		- Perform MCMC sampling with 9-mer fragments from kT=100 to kT=1 (assembly stage)
		- Perform MCMC sampling with 3-mer fragments from kT=1 to kT=0.1 (refinement stage)
		- Take best (lowest-energy) structure after refinement and perform energy minimization (see utils.relax)
		- Log energy and RMSD to native structure after minimization
	3. Visualize lowest-RMSD structure in PyMol

"""

from pyrosetta import *

init(extra_options='-mute all -constant_seed')
from FragmentSampler import MCMCSampler

import argparse
import os
from typing import Union


def formatArgs(arg: Union[list, str, int], convertfn: callable = int):
    if isinstance(arg, list):
        return convertfn(arg[0])
    elif isinstance(arg, str):
        return convertfn(arg)
    else:
        return arg


def main():
    # parsing arguments
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument("--fasta", metavar='--fasta', type=str, nargs=1, dest="fasta",
                        default=None, help=".fasta file containing sequence")
    parser.add_argument("--logdir", metavar="--logdir", type=str, nargs=1, dest="logdir",
                        default='.', help="directory to save all log files")
    parser.add_argument("--nsims", metavar="--nsims", type=int, nargs=1, dest="nsims",
                        default=1, help="number of simulations")
    parser.add_argument("--nfrags", metavar="--nfrags", type=int, nargs=1, dest="nfrags",
                        default=3, help="number of fragments to sample from at each iteration")
    parser.add_argument("--anneal_rate", metavar="--anneal_rate", type=float, nargs=1, dest="anneal_rate",
                        default=0.999, help="temperature annealing parameter")
    # getting system arguments
    args = parser.parse_args()
    nsims = formatArgs(args.nsims, int)
    # format N, could be list or string
    nfrags = formatArgs(args.nfrags, int)
    anneal_rate = formatArgs(args.anneal_rate, float)
    # start simulation
    seed = 9000
    result = [("sim_number", "energy", "rmsd")]
    for sim in range(1, nsims + 1):
        mcmc = MCMCSampler(args.fasta[0], args.logdir[0], sample_size=nfrags, annealing_rate=anneal_rate)
        result.append(mcmc.simulate(sim_index = sim, seed=seed + sim))
    with open(os.path.join(args.logdir[0], "simulation_summary.txt"), "w") as f:
        f.write('\n'.join('\t'.join(str(x) for x in i) for i in result))


if __name__ == '__main__':
    main()
