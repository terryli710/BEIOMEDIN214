"""
This file contains the main fragment sampling class, which performs a Monte Carlo simulated annealing procedure to fold a protein.
"""

from pyrosetta import *

init(extra_options='-mute all  -constant_seed')
from Bio.SeqIO import parse
import math
import utils
from Protein import Protein
from FragmentSet import FragmentSet
from typing import Union, Tuple
import os
import random


class MCMCSampler(object):

    def __init__(self, fasta: str, logdir: str = None, start_pdb: Union[str, None] = None, sample_size=50,
                 annealing_rate=0.999):
        """
        Initializing a MCMC sampler for certain protein
        The score function is given to you (Rosetta centroid score function)
        :param fasta: name fasta file
        :param sample_size: size of candidate fragments in each position
        attributes:
            scorefxn: score function, callable
            target_pose: goal pose, rosetta pose object
            protein: Protein class
            fragment_set: FragmentSet class
            self.candidate_frag: store calculated candidate fragments
            self.mers: k-mer
            self.temp: current temperature
            self.t_end: set ending temperature
        """
        ## 0
        # set log
        self.logdir = logdir
        ## 1
        # set score function
        self.scorefxn = create_score_function('score3')
        ## 2
        # read pdb file (goal position)
        # pose_from_pdb doesn't take absolute dir
        self.protein_name = fasta.split('.')[0]
        self.target_pose = pose_from_pdb(self.protein_name + '.pdb')
        ## 3
        # read fasta file (protein)
        fasta_path = os.path.join(self.protein_name + '.fasta')
        iter = parse(fasta_path, 'fasta')
        seq = next(iter)
        # initialize protein, either from seq or from start_pdb
        if start_pdb:
            os.chdir(os.path.dirname(start_pdb))
            self.protein = Protein(pose=pose_from_pdb(os.path.basename(start_pdb)))
        else:
            self.protein = Protein(sequence=seq.seq._data)
        # store initial pdb
        self.initial_protein = Protein(sequence=seq.seq._data)
        ## 4
        # get fragment set
        self.fragment_set = {"9mers": FragmentSet(os.path.join(self.protein_name + "_9mers.frag"),
                                                  os.path.join(self.protein_name + "_9mers.rmsd")),
                             "3mers": FragmentSet(os.path.join(self.protein_name + "_3mers.frag"),
                                                  os.path.join(self.protein_name + "_3mers.rmsd"))}
        ## 5
        # parametrize candidate_frag_dict
        self.candidate_frag = {"9mers": {}, "3mers": {}}
        for pos in range(1, self.fragment_set["9mers"].length + 1):
            self.candidate_frag["9mers"][pos] = self.fragment_set["9mers"].get_lowRMS_fragments(pos, sample_size)
        for pos in range(1, self.fragment_set["3mers"].length + 1):
            self.candidate_frag["3mers"][pos] = self.fragment_set["3mers"].get_lowRMS_fragments(pos, sample_size)
        ## 6
        # set temperature
        self.temp = 100
        self.t_end = 0.1
        ## 7
        # set anneal rate
        self.annealing_rate = annealing_rate
        return

    def compute_energy(self, protein: Union[Protein, None] = None) -> float:
        """
        compute energy of protein.
        Hint: look at utils.py
        --------
        Params:
            - protein (Protein object): protein to score
        Return:
            - energy of conformation (float)
        """
        # NOTE: score_pose cannot take absolute directory
        if protein:
            return utils.score_pose(protein.pose, self.scorefxn)
        else:
            return utils.score_pose(self.protein.pose, self.scorefxn)

    def perturb_fragment(self, pos: int, mer: str = "9mers",
                         protein: Union[Protein, None] = None) -> Tuple[
        Protein, int]:  # you may want to add more arguments
        """
        Sample from possible fragments for a position, and replace torsion angles of that fragment in the protein.
        Store fragment candidate at certain position (call get_lowRMS just once.)
        :param protein: optional parameter, if none, use self.protein
        :param pos: position to change
        :param mer: mode of function, either "3mers" or "9mers"
        :return: new Protein with updated angles
        """
        # set a new_pose (protein)
        if not protein:
            new_protein = Protein(pose=self.protein.pose)
        else:
            new_protein = Protein(pose=protein.pose)
        # sample candidate fragment
        random_index = random.randint(0, len(self.candidate_frag[mer][pos]) - 1)
        frag_chosen = self.candidate_frag[mer][pos][random_index]
        frag_index = self.fragment_set[mer].findFragIndex(pos, frag_chosen)
        # insert this fragment and return
        if mer == "9mers":
            frag_length = 9
        else:
            frag_length = 3
        for i in range(frag_length):
            new_protein.set_torsion(pos + i, frag_chosen[i][0], frag_chosen[i][1])
        return new_protein, frag_index

    def metropolis_accept(self, new_protein: Protein) -> float:  # you may want to add more arguments
        """
        Calculate probability of accepting or rejecting move based on Metropolis criterion.
        :param new_protein: candidate protein to be calculated and compared
        :return: probability of accepting
        """
        delta_e = self.compute_energy(new_protein) - self.compute_energy()
        # formula: if delta_E > 0: exp(-delta_E/kT)
        return math.exp(-delta_e / self.temp) if delta_e > 0 else 1

    def anneal_temp(self) -> bool:
        """
        Anneal temperature using exponential annealing schedule.
        Consider kT to be a single variable (i.e. ignore Boltzmann constant)
        :return whether it reached the threshold
        """
        assert self.temp > self.t_end, "Temperature has reached threshold"
        self.temp *= self.annealing_rate
        if self.temp <= self.t_end:
            return True
        else:
            return False

    def step(self, verbose=0) -> bool:
        """
        Take a single MCMC step. Each step should do the following:
        1. sample position in chain
            - Note: think about positions you can sample a k-mer fragment from. 
              For example, you cannot sample from position 1 because there is no phi angle
        2. sample fragment at that position and replace torsions in a *copied version* of the protein
        3. measure energy after replacing fragment
        4. accept or reject based on Metropolis criterion
            - if accept: incorporate proposed insertion and anneal temperature
            - if reject: sample new fragment (go to step 3)
        """
        accept = 0
        i = 0
        done = False
        if self.temp > 1:
            mer_str = "9mers"
        else:
            mer_str = "3mers"
        # 1. sample position in chain (e.g. len=10, 3-mers, should sample {1,...,7})
        sampled_pos = random.randint(1, self.fragment_set[mer_str].length)
        # get number of frag in this position
        pool_size = len(self.candidate_frag[mer_str][sampled_pos])
        sampled_set = set()
        # if accepted or sampled all frags and cannot decide, keep going
        while not accept and len(sampled_set) < pool_size:
            # 2. replace torsions in a *copied version* of the protein
            new_protein, index = self.perturb_fragment(sampled_pos, mer=mer_str)
            # add to set
            sampled_set.add(index)
            # 3. 4. measure energy and decide
            prob = self.metropolis_accept(new_protein)
            accept = random.uniform(0, 1) < prob
            if accept:
                # incorporate proposed insertion and anneal temperature
                self.protein = new_protein
                done = self.anneal_temp()
            # if reject: sample new fragment (go to step 2)
            i += 1
        if verbose:
            if accept:
                print("sampled position = {}, take {} iter to finish, prob is {}".format(sampled_pos, i, prob))
            elif len(sampled_set) == pool_size:
                print("sampled position = {}, didn't accept any frags".format(sampled_pos))
        return done

    def savelog(self, log: dict, file_name: str) -> None:
        """
        save log of sim
        :param log: log information
        :param file_name: saved path
        """
        saved_log = "iteration" + "\t" + "\t".join(log.keys()) + "\n"
        iter = 1
        for row in range(len(log["energy"])):
            saved_log += str(iter) + "\t"
            saved_log += "\t".join(str(log[key][iter - 1]) for key in log.keys())
            saved_log += "\n"
            iter += 1
        with open(file_name, "w") as f:
            f.write(saved_log)

    def storeSim(self, best_pdb: Protein, log: dict, sim_index: int) -> Tuple[str, str, int]:
        """
        Store best pdb and log text file to log
        :param best_pdb: the structure to store as "best.pdb"
        :param log: log dict to store
        :param sim_index: int of simulation number
        :return path to sim folder, path to log folder, sim
        """
        # dealing with paths
        if not self.logdir:
            cur_dir = os.getcwd()
            log_folder_name = self.protein_name + "_log"
            log_folder_path = os.path.join(cur_dir, log_folder_name)
        else:
            log_folder_path = self.logdir
        if not os.path.exists(log_folder_path):
            os.mkdir(log_folder_path)
        sim_folder_name = "sim_" + self.__toStr__(sim_index)
        sim_folder_path = os.path.join(log_folder_path, sim_folder_name)
        # avoid path exist error
        if not os.path.exists(sim_folder_path):
            os.mkdir(sim_folder_path)
        # store things
        # 1. initial pdb
        self.initial_protein.save_pdb(os.path.join(sim_folder_path, "initial.pdb"))
        # 2. target pdb
        target_protein = Protein(pose=self.target_pose)
        target_protein.save_pdb(os.path.join(sim_folder_path, "target.pdb"))
        # 3. best pdb
        best_pdb.save_pdb(os.path.join(sim_folder_path, "best.pdb"))
        # 6. log.txt
        self.savelog(log, os.path.join(sim_folder_path, sim_folder_name + "_log.txt"))
        return sim_folder_path, log_folder_path, sim_index

    @staticmethod
    def __toStr__(integer) -> str:
        """
        convert integer to formatted string
        :param integer: integer to be converted
        :return: string
        """
        if integer < 10:
            return "0" + str(integer)
        else:
            return str(integer)

    def simulate(self, sim_index: int, seed: int = 9001) -> Tuple[int, float, float]:
        """
        Run full MCMC simulation from start_temp to end_temp.
        Be sure to save the best (lowest-energy) structure, so you can access it after.
        It is also a good idea to track certain variables during the simulation (temp, energy, and more).
        :param sim_index: simulation index
        :param seed: int
        :return: log information
        """
        random.seed(seed)
        log = {"temperature": [], "energy": []}
        while True:
            if self.step():
                # store best pdb
                sim_path, log_folder_path, sim_index = self.storeSim(self.protein, log, sim_index)
                # calculate relaxed, cd to the folder!!
                cur_dir = os.getcwd()
                os.chdir(sim_path)
                protein, rmsd, score = utils.relax("best.pdb",
                                                   "target.pdb")
                os.chdir(cur_dir)
                break
            else:
                # keep track of log
                log['temperature'].append(self.temp)
                log['energy'].append(self.compute_energy())
        return sim_index, score, rmsd
