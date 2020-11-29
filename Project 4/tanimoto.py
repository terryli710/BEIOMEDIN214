''' This file outputs a csv file that contains <drugidA>,<drugidB>,<Tanimoto_score>,<share_target>'''

import sys
from collections import defaultdict
from itertools import combinations
from typing import Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from chemoUtils import Targets, Drug


class Tanimoto(object):

    def __init__(self, drug_file, targets_file):
        '''
        initialize instance
        :param drug_file: drug_file path
        :param targets_file: targets file path
        '''
        self.drug = Drug(drug_file)
        self.targets = Targets(targets_file)
        self.result_dict = defaultdict(lambda: {})
        self.result_filled = False

    def getDrugPair(self) -> set:
        ''' return all drug pairs '''
        return set(combinations(self.drug.getAllDrug(), 2))

    def calcTanimoto(self, drugA: str, drugB: str) -> float:
        '''
        calculate tanimoto score using the equation
        :param drugA: db_id of drug A
        :param drugB: db_id of drug B
        :return: float of tanimoto score
        '''
        maccs_a, maccs_b = self.drug.getMaccs(drugA), self.drug.getMaccs(drugB)
        return len(maccs_a & maccs_b) / len(maccs_a | maccs_b)

    def targetIntersect(self) -> set:
        '''
        return return list of all drug pair that have share targets
        :return: list of shared targets
        '''
        result = []
        for protein in self.targets.getAllProtein():
            bind_drug_set = self.targets.getLigands(protein)
            result += list(combinations(bind_drug_set, 2))
        return set(result)

    def fillResult(self) -> None:
        ''' function that fills self.result_dict '''
        drug_pair_set = self.getDrugPair()
        shared_drug_pairs = self.targetIntersect()
        print("Filling results")
        for drug_pair in drug_pair_set:
            t_score = self.calcTanimoto(drug_pair[0], drug_pair[1])
            shared_targets = 1 if drug_pair in shared_drug_pairs else 0
            self.result_dict[drug_pair] = {'tanimoto_score': t_score,
                                           'shared_target': int(bool(shared_targets))}
        self.result_filled = True
        print("Results filled")
        return

    def getResultCsv(self) -> Union[None, pd.DataFrame]:
        ''' returns the pandas DataFrame for csv saving '''
        if not self.result_filled:
            print('Result not computed')
            return None
        else:
            drugA_lst, drugB_lst, t_score_lst, shared_target_lst = [], [], [], []
            for drug_pair in self.result_dict.keys():
                drugA_lst.append(drug_pair[0])
                drugB_lst.append(drug_pair[1])
                t_score_lst.append(self.result_dict[drug_pair]['tanimoto_score'])
                shared_target_lst.append(self.result_dict[drug_pair]['shared_target'])
            result_df = pd.DataFrame(data={'drug_a': drugA_lst,
                                           'drug_b': drugB_lst,
                                           'tanimoto_score': t_score_lst,
                                           'shared_target': shared_target_lst})
            return result_df

    def saveCsv(self, save_root) -> None:
        '''
        Save csv to save root
        :param save_root: path and csv name to save
        '''
        result_df = self.getResultCsv()
        if isinstance(result_df, pd.DataFrame):
            result_df.to_csv(save_root, index=False, header=False, float_format='%.6f')
            print('Saved to ', save_root)

    def drawHistograms(self):
        ''' Draw histograms by calling __drawHist__ '''
        sunitID = 'yyhhli'
        # plot all
        result_df = self.getResultCsv()
        self.__drawHist__(result_df['tanimoto_score'], sunitID + ' All', 'all_tanimoto.png')
        # plot shared
        shared_df = result_df[result_df['shared_target'] == 1]['tanimoto_score']
        self.__drawHist__(shared_df, sunitID + ' Shared', 'shared_tanimoto.png')
        # plot notshared
        not_shared_df = result_df[result_df['shared_target'] == 0]['tanimoto_score']
        self.__drawHist__(not_shared_df, sunitID + ' Not Shared', 'notshared_tanimoto.png')

    @classmethod
    def __drawHist__(cls, df_t_score: pd.DataFrame, title: str, save_name: str) -> None:
        '''
        Draw histogram
        :param df_t_score: tanimoto score as a coloumn of pandas dataframe
        :param title: title of the plot
        :param save_name: saving name of the plot
        '''
        plt.hist(df_t_score, bins=cls.calcNBins(df_t_score))
        plt.title(title)
        plt.xlabel('Tanimoto Score')
        plt.ylabel('Count')
        plt.savefig(save_name, dpi=400)
        plt.close()

    @classmethod
    def calcNBins(cls, df_t_score):
        ''' Calculating number of bins for histogram using scott'''
        return int(
            (np.max(df_t_score) - np.min(df_t_score)) / (3.5 * np.std(df_t_score) / (df_t_score.shape[0] ** (1 / 3))))


def main():
    ''' Main function, save csv and save plots'''
    drug_file, target_file, save_root = sys.argv[1:]
    t = Tanimoto(drug_file, target_file)
    t.fillResult()
    t.saveCsv(save_root)
    t.drawHistograms()


if __name__ == '__main__':
    main()
