"""
 This file implements GSEA (as discussed in lecture) to identify gene sets that are significantly enriched
 in pathological samples to get an idea of potentially activated pathways involved in the disease.
 We will be using pathways from the KEGG database.
"""

import os
import sys
import pandas as pd
import numpy as np
from typing import List, Union


class GSEA:
    def __init__(self) -> None:
        """
        Initializing GSEA instance
        """
        self.expfile = None
        self.sampfile = None
        self.genesets = None

    def load_data(self, expfile: str, sampfile: str, genesets: str) -> None:
        """
        take the file paths to the expression, sample and gene set data,
        read them in and store within the GSEA instance
        :param expfile: expression file path
        :param sampfile: sample file path
        :param genesets: gene set data file path
        """
        # if either nothing inputs, or can't find file, pretend nothing happened
        if not os.path.isfile(expfile) or \
                not os.path.isfile(sampfile) or \
                not os.path.isfile(genesets):
            print("File not found")
            return
        self.expfile = pd.read_csv(expfile, delimiter='\t', header=0, index_col=0)
        self.sampfile = pd.read_csv(sampfile, delimiter='\t', header=None)
        with open(genesets, 'r') as file:
            content = file.read()
        # format geneset as list of list, throw empty lines
        self.genesets = [line.split('\t') for line in content.split('\n') if line != '']
        return

    def get_gene_rank_order(self, permute: bool = False) -> List[str]:
        """
        return a list of all genes (as strings) ranked by their logFC between patient and control,
        with the gene with the highest logFC ordered first.
        :param permute: boolean value, "true" used in calculating p-value
        :return: list of all gene ranked by logFC = log(P) - log(C)
        """
        # get index of positive and negative cases (optional permutation)
        if not permute:
            pos_index = list(self.sampfile[self.sampfile[1] == 1][0])
            neg_index = list(self.sampfile[self.sampfile[1] == 0][0])
        else:
            index = list(np.random.permutation(self.sampfile[0]))
            num_pos = np.sum(self.sampfile[1] == 1)
            pos_index = index[:int(0.5 * num_pos)]
            neg_index = index[int(0.5 * num_pos):]
        # get mean of pos and neg classes
        pos_mean = np.mean(self.expfile.loc[:, pos_index], axis=1)
        neg_mean = np.mean(self.expfile.loc[:, neg_index], axis=1)
        # rank the gene (NOTE: it should be largest logFC first, so sort by negative value)
        return list(self.expfile.index[(-(pos_mean - neg_mean)).argsort()])

    def get_enrichment_score(self, geneset: str, permute=False) -> float:
        """
        return the enrichment score, a float correct to two decimal places, for a given gene set,
        such as `KEGG_CITRATE_CYCLE_TCA_CYCLE`.
        This method should run get_gene_rank_order at some point to initiate enrichment calculations.
        :param permute: boolean value, "true" used in calculating p-value
        :param geneset: name of gene set of interest
        :return: the enrichment score, a float correct to two decimal places
        """
        # ranked gene list
        gene_list = self.get_gene_rank_order(permute=permute)
        # 1. calculate up and down score (step size)
        # get the genes in gene set (filtered to be in our files)
        genes = [gene for gene in self.__getGeneSet__(geneset) if gene in self.expfile.index]
        # num of total genes
        nt = self.expfile.shape[0]
        # num of genes in gene set
        ng = len(genes)
        up_score = np.sqrt((nt - ng) / ng)
        down_score = - np.sqrt(ng / (nt - ng))
        # 2. moving down the ranked list of genes.
        bb_score = [0]
        current_score = 0
        for gene in gene_list:
            # if target gene, add score
            if gene in genes:
                current_score += up_score
            # if no, subtract score
            else:
                current_score += down_score
            bb_score.append(current_score)
        # 3. find supremum score (not taking absolute) (round to two decimal place)
        return round(max(bb_score), 2)

    def __getGeneSet__(self, geneset: str) -> Union[List[str], None]:
        """
        return the genes in the specified gene set name
        :param geneset: name of gene set of interest
        :return: list of names of genes
        """
        try:
            gene_set_names = [gs[0] for gs in self.genesets]
            index = gene_set_names.index(geneset)
            return self.genesets[index][2:]
        except ValueError:
            print('Can\'t find specified gene set.')
            return

    def getSignificance(self, geneset: str, permutation_num: int = 100) -> float:
        """
        Calculate p-value for a specific geneset
        :param geneset: name of gene set of interest
        :return: p-value of significance
        """
        score = self.get_enrichment_score(geneset)
        random_score = []
        for i in range(permutation_num):
            random_score.append(self.get_enrichment_score(geneset, permute=True))
        pvalue = np.mean(np.array(random_score) > score)
        return pvalue

    def get_sig_sets(self, p: float, permutation_num: int = 100) -> List[str]:
        """
         return the list of significant gene sets (as strings),
         at a corrected threshold of p, by name.
        :param p: corrected threshold of p
        :return: list of significant gene sets
        """
        geneset_names = [i[0] for i in self.genesets]
        pvalue = []
        for geneset_name in geneset_names:
            # print("working on set ", geneset_name)
            p_for_set = self.getSignificance(geneset_name, permutation_num=permutation_num)
            # print("for set {}, p is {}".format(geneset_name, p_for_set))
            pvalue.append(p_for_set)
        return [geneset_names[i] for i in range(len(geneset_names)) if pvalue[i] < p]


def main(expfile, sampfile, keggfile) -> None:
    """
    Main function for command line usage (should be modified before using)
    :param expfile: expression file path
    :param sampfile: sample file path
    :param genesets: gene set data file path
    :return: None
    """
    gsea = GSEA()
    gsea.load_data(expfile, sampfile, keggfile)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3])
    # gsea = GSEA()
    # gsea.load_data('GSE25628_filtered_expression.txt', 'GSE25628_samples.txt', 'c2.cp.kegg.v6.2.symbols.filtered.gmt')
    # gsea.get_sig_sets(p=0.05)
