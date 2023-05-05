import networkx as nx
import pandas as pd
import numpy as np
import scipy
import matplotlib.pyplot as plt
import gseapy as gp
import seaborn as sns

import itertools
import sys
import pprint

def compute_corr(
    gene_set_paths,
    is_mcp=False,
    progress_increment=1000):
    """PLAN
    1. Auxin and Ethylene GO next to each other as super node anchors
    2. The nodes are the individual XM_
    3. The edges are their correlations to each other. Have to do threshold
    4. (maybe) remove Nodes that don't have any edges
    5. Identify the nodes that are in both Auxin and Ethylene

    """
    df_list = [pd.read_excel(f"./data/processed/03_gene_sets/{gene_set}.xlsx") for gene_set in gene_set_paths.values()]
    df_combined = pd.concat([*df_list], axis=0)
    num_entries = np.size(df_combined,0)
    corr_matrix = np.zeros((num_entries*num_entries))

    expressions = df_combined.iloc[:,10:15] if is_mcp else df_combined.iloc[:,4:10]

    for i, (gene_a, gene_b) in enumerate(itertools.combinations(expressions.values.tolist(), 2)):
        scipy.stats.pearsonr(gene_a, gene_b)
        corr_matrix[i] = scipy.stats.pearsonr(gene_a, gene_b).statistic
        if i % progress_increment == 0:
            print(f"Progress: ({i} of {num_entries**2})", end="\r")
            sys.stdout.write("\033[K")

    corr_matrix = np.reshape(corr_matrix, (num_entries, num_entries))
    output_name = "mcp" if is_mcp else "ctrl"
    np.save(f"./data/processed/03_gene_sets/corr_matrix_{output_name}.npy", corr_matrix)

def graph_corr(corr_matrix_file, title):
    corr_matrix = np.load(corr_matrix_file)
    auxin = pd.read_excel("./data/processed/03_gene_sets/P:auxin-activated signaling pathway.xlsx")
    ethylene = pd.read_excel("./data/processed/03_gene_sets/P:ethylene-activated signaling pathway.xlsx")
    df_combined = pd.concat([auxin, ethylene], axis=0, ignore_index=True)
    

    filter_matrix = abs(corr_matrix) > .999
    # np_combined = df_combined.to_numpy()

    G = nx.Graph()

    # Plot nodes
    print(len(df_combined))
    G.add_node("Auxin", go_id={"auxin"}, go_node=True)
    
    G.add_node("Ethylene", go_id={"ethylene"}, go_node=True)
    eth_nodes = [None] * np.size(ethylene,0)
    parent_edges = []
    for i, row in enumerate(ethylene["accession_number"]):
        # if df_combined.loc[i,"accession_number"] in set(ethylene["accession_number"]):
        G.add_node(row, go_id={"ethylene"}, go_node=False)
        G.add_edge(row, "Ethylene", weight=2)
        parent_edges.append((row, "Ethylene"))

    for i, row in enumerate(auxin["accession_number"]):
        # if df_combined.loc[i,"accession_number"] in set(auxin["accession_number"]):
        if G.has_node(row):
            G.nodes[row]["go_id"].add("auxin")
        else:
            G.add_node(row, go_id={"auxin"}, go_node=False)
        G.add_edge(row, "Auxin", weight=2)

        parent_edges.append((row, "Auxin"))

    edgelist = []
    
    for i, row in enumerate(corr_matrix):
        for j, col in enumerate(row):
            if (filter_matrix[i,j]) and (i != j):
                u_node = df_combined.loc[i, "accession_number"]
                v_node = df_combined.loc[j, "accession_number"]
                G.add_edge(u_node, v_node)
                edgelist.append((u_node, v_node))

        if i % 100 == 0:
            print(f"We are on row {i} of {np.size(corr_matrix,0)}")

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
        )
    
    for go in ["auxin", "ethylene"]:
        nx.draw_networkx_nodes(
            G, pos, 
            nodelist=[node for node, attr in G.nodes.data() if go in attr["go_id"] and not attr["go_node"]], 
            node_color=("red" if go=="ethylene" else "blue"), 
            node_size=10, 
            alpha=0.25)
    
    nx.draw_networkx_edges(G, pos, edgelist=edgelist, edge_color="purple", alpha=0.25)
    

    print(f"Num edges: {G.number_of_edges()}")
    print(f"Num nodes: {G.number_of_nodes()}")
    plt.title(title)
    plt.savefig(f"./figures/{title}", dpi=1000)
    plt.show()

if __name__ == "__main__":
    # compute_corr()
    graph_corr(
    "./data/processed/03_gene_sets/corr_matrix_ctrl.npy", 
    title="Ethylene-Auxin Crosstalk Network After Ethylene-Receptor Inhibition")