''' This file uses networkx to generate a visualization  network using 'network_edgelist.txt' '''

import sys

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd


def main() -> None:
    ''' Main function for plotting network '''
    network_edgelist_file, protein_nodes_file, output_path = sys.argv[1:]
    # load files
    network_edgelist = pd.read_csv(network_edgelist_file, sep=' ', header=None)
    network_edgelist.columns = ['node1', 'node2']
    protein_nodes = pd.read_csv(protein_nodes_file)
    # filter protein_nodes to only include nodes that are in network_edgelist
    all_nodes = list(network_edgelist['node1']) + list(network_edgelist['node2'])
    protein_nodes = protein_nodes[protein_nodes['uniprot_accession'].isin(all_nodes)]
    # create node labels dictionary
    node_label_dict = {}
    for i, row in protein_nodes.iterrows():
        node_label_dict[row['uniprot_accession']] = row['uniprot_id']
    # create color mapping
    color_map = {"bp": "red",
                 "bp;cholesterol": "green",
                 "bp;cholesterol;diabetes": "blue",
                 "bp;diabetes": "purple"}
    # create graph
    G = nx.from_pandas_edgelist(network_edgelist, 'node1', 'node2')
    # draw graph
    plt.figure(figsize=(8, 8))
    pos = nx.spring_layout(G)
    options = {"node_size": 500, "alpha": 0.8}
    for label, color in color_map.items():
        node_lst = list(protein_nodes.loc[protein_nodes.indications == label, 'uniprot_accession'])
        nx.draw_networkx_nodes(G, pos, nodelist=node_lst, node_color=color, **options)
    nx.draw_networkx_edges(G, pos)
    nx.draw_networkx_labels(G, pos, labels=node_label_dict)
    plt.tight_layout()
    plt.savefig('network.png', format='PNG', dpi=150)


if __name__ == '__main__':
    main()
