"""
This file contains the main fragment sampling class, which performs a Monte Carlo simulated annealing procedure to fold a protein.
"""

from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.scoring import *
init(extra_options='-mute all  -constant_seed')
from Bio.SeqIO import parse
import numpy as np
import utils
from Protein import Protein
from FragmentSet import FragmentSet

import os
import random




class MCMCSampler(object):
    def __init__(self, protein_name: str, data_dir: str, mers: int):
        """
        Initializing a MCMC sampler for certain protein
        The score function is given to you (Rosetta centroid score function)
        :param protein_name: name of that protein
        """
        # set score function
        self.scorefxn = create_score_function('score3')
        # read pdb file (goal position)
        self.target_pose = pose_from_pdb(os.path.join(data_dir, protein_name+'.pdb'))
        # read fasta file (protein)
        fasta_path = os.path.join(data_dir, protein_name + '.pdb')
        iter = parse(fasta_path)
        seq = next(iter)
        self.protein = Protein(sequence=seq.seq._data)
        self.fragment_set = FragmentSet(os.path.join(data_dir, protein_name + "{}mers.frag".format(str(mers))),
                                        os.path.join(data_dir, protein_name + "{}mers.rmsd".format(str(mers))))
        return


    def compute_energy(self, protein: Protein) -> float:
        """
        compute energy of protein.
        Hint: look at utils.py
        --------
        Params:
            - protein (Protein object): protein to score
        Return:
            - energy of conformation (float)
        """        
        return utils.score_pose(protein, self.scorefxn)

    def perturb_fragment(self, protein: Protein, pos: int) -> Protein: # you may want to add more arguments
        """
        TODO: Sample from possible fragments for a position, and replace torsion angles of that fragment in the protein.
        ---------
        Params:
            - TODO
        Returns:
            - TODO
        """
        pass


    def metropolis_accept(self): # you may want to add more arguments
        """
        TODO: Calculate probability of accepting or rejecting move based on Metropolis criterion.
        --------
        Params:
            - TO DO
        Returns:
            - TO DO
        """
        pass

    def anneal_temp(self):
        """
        TODO: Anneal temperature using exponential annealing schedule. Consider kT to be a single variable (i.e. ignore Boltzmann constant)
        --------
        Params:
            - TO DO
        Returns:
            - TO DO
        """
        pass

    def step(self):
        """
        TODO: Take a single MCMC step. Each step should do the following:
        1. sample position in chain
            - Note: think about positions you can sample a k-mer fragment from. 
              For example, you cannot sample from position 1 because there is no phi angle
        2. sample fragment at that position and replace torsions in a *copied version* of the protein
        3. measure energy after replacing fragment
        4. accept or reject based on Metropolis criterion
            - if accept: incorporate proposed insertion and anneal temperature
            - if reject: sample new fragment (go to step 3)
        """
        pass


    def simulate(self):
        """
        TODO: Run full MCMC simulation from start_temp to end_temp.
        Be sure to save the best (lowest-energy) structure, so you can access it after.
        It is also a good idea to track certain variables during the simulation (temp, energy, and more).
        -------- 
        Params:
            - TO DO
        Returns:
            - TO DO
        """
        pass


def debug():
    from pathlib import Path
    cur_path = os.getcwd()
    project_dir = Path(cur_path).parent
    data_dir = os.path.join(cur_path, "starter_data")
    doc_name = "1fw4.fasta"
    iterator = parse(os.path.join(data_dir, doc_name), "fasta")
    seq = next(iterator)
    p = Protein(sequence=seq.seq._data)
    frag = "helix_9mers.frag"
    rmsd = "helix_9mers.rmsd"
    fs = FragmentSet(os.path.join(data_dir, frag), os.path.join(data_dir, rmsd))
    fs.get_lowRMS_fragments(3, 3)
    mcmc = MCMCSampler("helix", data_dir, 9)
    mcmc.perturb_fragment(p, pos=3)

if __name__ == '__main__':
    debug()


    



