''' This file generates node list as a txt file for network plotting '''

import sys
from itertools import combinations

import pandas as pd

from pvalue import Similarity


def main() -> None:
    ''' Main function to generate 'network_edgelist.txt' '''
    # get arguments
    drug_file, targets_file, protein_node_file = sys.argv[1:]
    # set up instance
    sim = Similarity(drug_file, targets_file)
    # load file
    protein_nodes = pd.read_csv(protein_node_file)
    # generate all pairs
    protein_pairs_lst = list(combinations(protein_nodes['uniprot_accession'], 2))
    # calculate all p-values
    result = []
    print("Start calculating p-values")
    for protein_pair in protein_pairs_lst:
        p_value = sim.calcPValue(protein_pair[0], protein_pair[1], num=500, random_seed=214)
        # if significant
        if p_value <= 0.05:
            result.append(protein_pair)
    # save to file
    result_str = '\n'.join([' '.join(item) for item in result])
    with open('network_edgelist.txt', 'w') as f:
        f.write(result_str)
    print("Saved to network_edgelist.txt")
    return


if __name__ == '__main__':
    main()