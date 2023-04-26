"""This module handles in-house gene set along with enrichment analysis

By creating gene sets based on the the GO:0009755 (hormone-mediated signaling pathway),
we can identify if there are any entities that are part of multiple hormone-mediated 
pathways. These are our "cross-talk regions" of interest (de Anda-Jáuregui et al. 
BMC Systems Biology).

"""

import networkx as nx
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt

import itertools
import sys

def compute_corr():
    """PLAN

    1. Auxin and Ethylene GO next to each other as super node anchors
    2. The nodes are the individual XM_
    3. The edges are their correlations to each other. Have to do threshold
    4. (maybe) remove Nodes that don't have any edges
    5. Identify the nodes that are in both Auxin and Ethylene

    """
    auxin = pd.read_excel("./data/processed/03_gene_sets/P:auxin-activated signaling pathway.xlsx")
    ethylene = pd.read_excel("./data/processed/03_gene_sets/P:ethylene-activated signaling pathway.xlsx")
    df_combined = pd.concat([auxin, ethylene], axis=0)
    num_entries = np.size(df_combined,0)
    corr_matrix = np.zeros((num_entries*num_entries))

    # expressions = df_combined.iloc[:,4:10] # JUST CONTROLS
    expressions = df_combined.iloc[:,10:15] # JUST MCP
    # print(expressions.iloc[1,:])

    for i, (gene_a, gene_b) in enumerate(itertools.combinations(expressions.values.tolist(), 2)):

        scipy.stats.pearsonr(gene_a, gene_b)
        corr_matrix[i] = scipy.stats.pearsonr(gene_a, gene_b).statistic

        if i % 1000 == 0:
            print(f"Progress: ({i} of {num_entries**2})", end="\r")
            sys.stdout.write("\033[K")

    corr_matrix = np.reshape(corr_matrix, (num_entries, num_entries))
    np.save("./data/processed/03_gene_sets/corr_matrix.npy", corr_matrix)


def graph_corr(nodes, corr_matrix_file):
    corr_matrix = np.load(corr_matrix_file)
    auxin = pd.read_excel("./data/processed/03_gene_sets/P:auxin-activated signaling pathway.xlsx")
    ethylene = pd.read_excel("./data/processed/03_gene_sets/P:ethylene-activated signaling pathway.xlsx")
    df_combined = pd.concat([auxin, ethylene], axis=0, ignore_index=True)
    

    filter_matrix = abs(corr_matrix) > .999
    # np_combined = df_combined.to_numpy()

    G = nx.Graph()

    # Plot nodes
    print(len(df_combined))
    G.add_node("Auxin", go_id={"auxin"})
    
    G.add_node("Ethylene", go_id={"ethylene"})
    eth_nodes = [None] * np.size(ethylene,0)
    parent_edges = []
    for i, row in enumerate(ethylene["accession_number"]):
        # if df_combined.loc[i,"accession_number"] in set(ethylene["accession_number"]):
        G.add_node(row, go_id={"ethylene"})
        G.add_edge(row, "Ethylene")
        parent_edges.append((row, "Ethylene"))

    for i, row in enumerate(auxin["accession_number"]):
        # if df_combined.loc[i,"accession_number"] in set(auxin["accession_number"]):
        if G.has_node(row):
            G.nodes[row]["go_id"].add("auxin")
        else:
            G.add_node(row, go_id={"auxin"})
        G.add_edge(row, "Auxin")

        parent_edges.append((row, "Auxin"))



    edgelist = []
    for i, row in enumerate(corr_matrix):
        for j, col in enumerate(row):
            if (filter_matrix[i,j]) and (i != j):
                u_node = df_combined.loc[i, "accession_number"]
                v_node = df_combined.loc[j, "accession_number"]
                # print(u_node)
                # print(type(v_node))
                G.add_edge(u_node, v_node)
                edgelist.append((u_node, v_node))

        if i % 100 == 0:
            print(f"We are on row {i} of {np.size(corr_matrix,0)}")

    # print(edgelist)
    G.remove_nodes_from(list(nx.isolates(G)))

    pos = nx.spring_layout(
        G, 
        pos={
            "Auxin": (20,0),
            "Ethylene": (0,0)
        },
        fixed=["Auxin", "Ethylene"],
        k=1,
        threshold=1e-8,
        weight=1
        )
    
    nx.draw_networkx_nodes(G, pos, nodelist=["Ethylene"], node_color="red", node_size=50, label="Ethylene")
    nx.draw_networkx_nodes(G, pos, nodelist=["Auxin"], node_color="blue", node_size=50, label="Auxin")
    for go in ["auxin", "ethylene"]:

        nx.draw_networkx_nodes(
            G, pos, 
            nodelist=[node for node, attr in G.nodes.data() if go in attr["go_id"]], 
            node_color=("red" if go=="ethylene" else "blue"), 
            node_size=10, 
            alpha=0.25)
    
    # nx.draw
    nx.draw_networkx_edges(G, pos, edgelist=parent_edges, edge_color="gray", alpha=0.05)
    nx.draw_networkx_edges(G, pos, edgelist=edgelist, edge_color="purple", alpha=0.25)
    


    print(f"Num edges: {G.number_of_edges()}")
    print(f"Num nodes: {G.number_of_nodes()}")

    plt.show()


if __name__ == "__main__":
    # compute_corr()
    graph_corr([], "./data/processed/03_gene_sets/corr_matrix.npy")
    # G = nx.Graph()
    # G.add_node("hi", age=3)
    # G.add_node("hi", name="kevin")
    # print(G.nodes["hi"])
    # a = np.array([
    #     [1, 2, 3],
    #     [5, 6, 7],
    # ])
    # for i in a:
    #     print(i)