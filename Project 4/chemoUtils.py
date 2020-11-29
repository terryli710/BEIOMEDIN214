""" This file provides util functions that are useful across code files to prevent duplicate codes"""

import os
from collections import defaultdict
from typing import List

import pandas as pd


class Drug(object):
    '''
    class for drug file loading and related operations
    '''

    def __init__(self, drug_file_path: str) -> None:
        '''
        Load file data
        :param drug_file_path: path of 'drug.csv'
        '''
        # loading file
        drug_file = pd.read_csv(os.path.join(os.getcwd(), drug_file_path))
        # creating data as a dictionary
        self.drug_dict = defaultdict()
        for i, row in drug_file.iterrows():
            maccs = set(int(i) for i in row['maccs'].split(' '))
            self.drug_dict[row['db_id']] = {'generic_name': row['generic_name'],
                                            'maccs': maccs}
        pass

    def getMaccs(self, drug_id: str) -> set:
        ''' returns maccs given drug id'''
        return self.drug_dict[drug_id]['maccs']

    def getGenericName(self, drug_id: str) -> str:
        ''' returns generic_name given drug id'''
        return self.drug_dict[drug_id]['generic_name']

    def getAllDrug(self) -> List:
        ''' returns list of all drug in list'''
        return list(self.drug_dict.keys())


class Targets(object):
    '''
    class for loading target file and provide protein-based dictionary for searching
    '''

    def __init__(self, target_file_path: str) -> None:
        '''
        load file data
        :param target_file_path: path of 'targets.csv'
        '''
        # loading file
        target_file = pd.read_csv(os.path.join(os.getcwd(), target_file_path))
        # create data dictionary that link uniprot_accession to id and set of targets
        self.protein_dict = defaultdict(lambda: {})
        for i, row in target_file.iterrows():
            # if key is empty dict
            if not self.protein_dict[row['uniprot_accession']]:
                self.protein_dict[row['uniprot_accession']] = {'uniprot_id': row['uniprot_id'],
                                                               'db_id_set': {row['db_id']}}
            # if key exists
            else:
                self.protein_dict[row['uniprot_accession']]['db_id_set'].add(row['db_id'])
        pass

    def getId(self, accession_name: str) -> str:
        '''
        get uniprot_id from accession
        :param accession_name: accession input
        :return: id string
        '''
        return self.protein_dict[accession_name]['uniprot_id']

    def getLigands(self, accession_name: str) -> set:
        '''
        get a set of ligands (drugs) that the protein binds to
        :param accession_name: accession input
        :return: set of string containing names of drug as 'DB_id'
        '''
        return self.protein_dict[accession_name]['db_id_set']

    def getAllProtein(self) -> list:
        '''get all protein accession as a list'''
        return list(self.protein_dict.keys())
