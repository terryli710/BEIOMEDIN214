''' This file prints the p-value for two proteins to indicate their similarity'''

from chemoUtils import Drug, Targets
from tanimoto import Tanimoto

import argparse
from itertools import product
from typing import Union

import numpy as np


class Similarity(object):

    def __init__(self, drug_file, targets_file):
        ''' initialize the instance'''
        self.tanimoto = Tanimoto(drug_file, targets_file)

    def calcT(self, set_a: Union[set, list], set_b: Union[set, list]) -> float:
        '''
        calculate T_summary for two ligand set
        :param set_a: ligand set a
        :param set_b: ligand set b
        :return: float of t summary score
        '''
        result = 0
        drug_pair_lst = [sorted(item) for item in list(product(set_a, set_b))]
        for drug_pair in drug_pair_lst:
            t_score = self.tanimoto.calcTanimoto(drug_pair[0], drug_pair[1])
            result += (t_score if t_score > 0.5 else 0)
        return result

    def calcTsummary(self, proteinA: str, proteinB: str) -> float:
        '''
        calculate T summary score for two protein
        :param proteinA: accession of protein a
        :param proteinB: accession of protein b
        :return: float of t summary score
        '''
        # get all ligands (drugs that proteins bind to)
        ligands_a = self.__getLigandsSet__(proteinA)
        ligands_b = self.__getLigandsSet__(proteinB)
        return self.calcT(ligands_a, ligands_b)

    def __getLigandsSet__(self, protein) -> set:
        ''' get a set of ligands for a particular protein '''
        return self.tanimoto.targets.getLigands(protein)

    def bootstrapT(self, len_set_a: int, len_set_b: int, num: int, random_seed: int) -> list:
        '''
        calculate a list of bootstrap T scores
        :param len_set_a: length of ligand set a
        :param len_set_b: length of ligand set b
        :return: list of bootstrapped t summary score
        '''
        result = []
        for i in range(num):
            np.random.seed(int(random_seed+i))
            # sample from all the drugs (ligands) that we have
            set_a_boot = list(np.random.choice(list(self.tanimoto.drug.getAllDrug()), size=len_set_a, replace=True))
            set_b_boot = list(np.random.choice(list(self.tanimoto.drug.getAllDrug()), size=len_set_b, replace=True))
            result.append(self.calcT(set_a_boot, set_b_boot))
        return result

    def calcPValue(self, proteinA: str, proteinB: str, num: int, random_seed: int) -> float:
        '''
        Calculate p-value for two protein using bootstrap
        :param proteinA: accession of protein a
        :param proteinB: accession of protein b
        :param num: number of iterations for bootstrapping
        :param random_seed: random seed used for sampling
        :return: p-value rounded to 6 decimal place
        '''
        t_summary = self.calcTsummary(proteinA, proteinB)
        t_boot_lst = self.bootstrapT(len(self.__getLigandsSet__(proteinA)),
                                     len(self.__getLigandsSet__(proteinB)),
                                     num=num,
                                     random_seed=random_seed)
        return round(sum([int(t_boot >= t_summary) for t_boot in t_boot_lst]) / num, 6)


def main():
    # parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', dest='num_iter', type=int, default=500)
    parser.add_argument('-r', dest='random_seed', type=int, default=214)
    parser.add_argument('files', type=str, nargs=4)
    args = parser.parse_args()
    drug_file, targets_file, proteinA, proteinB = args.files
    # initiate some class objects
    sim = Similarity(drug_file, targets_file)
    print("{0:.6f}".format(sim.calcPValue(proteinA, proteinB, args.num_iter, args.random_seed)))


if __name__ == '__main__':
    main()