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
import gseapy as gp
import seaborn as sns
from IPython.display import display

import itertools
import sys
import pprint

def compute_corr(
    gene_set_paths=[
        "P:auxin-activated signaling pathway.xlsx", 
        "P:ethylene-activated signaling pathway.xlsx"],
    is_mcp=False,
    progress_increment=1000):
    """PLAN
    1. Auxin and Ethylene GO next to each other as super node anchors
    2. The nodes are the individual XM_
    3. The edges are their correlations to each other. Have to do threshold
    4. (maybe) remove Nodes that don't have any edges
    5. Identify the nodes that are in both Auxin and Ethylene

    """
    df_list = [pd.read_excel(f"./data/processed/03_gene_sets/{gene_set}") for gene_set in gene_set_paths]
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

def graph_corr(nodes, corr_matrix_file, title):
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

def gsea(
    phenotype: list,
    gene_set_paths: dict[str, str],
    is_ctrl: bool,
    title_gsea: str,
    title_heatmap: str,
    plot_enabled: bool=True):
    """Outputs GSEA data

    The steps involved with GSEA go as follows:
        1. Calculation of the Enrichment Score
        2. Estimation of the Sign
    """
    def retrieve_expression_data(
        gene_set_paths
        ) -> pd.DataFrame:
        dir_out = "./data/processed"
        df_list = []

        for gene_set, gene_path in gene_set_paths.items():
            df_gene_set = pd.read_excel(f"./data/processed/03_gene_sets/{gene_path}.xlsx")
            df_gene_set["gene_set"] = pd.Series(gene_set, index=range(np.size(df_gene_set,0)))
            df_list.append(df_gene_set)

        df_combined = pd.concat([*df_list], axis=0, ignore_index=True)
        df_combined = df_combined.drop_duplicates("accession_number")

        return df_combined
    
    def get_gene_sets(
        gene_set_paths: dict[str, str],
        path_base_dir="./data/processed/03_gene_sets/"
        ) -> dict[str, set]:
        """Returns a dictionary mapping "gene_set_name" to a Python set of genes

        This function prepares data for the input of gsea.calculate_enrichment_score
        """
    
        result = dict()
        for gene_set_name, gene_set_path in gene_set_paths.items():
            df_gene_set = pd.read_excel(f"{path_base_dir}{gene_set_path}.xlsx")
            result[gene_set_name] = set(df_gene_set["accession_number"])
        
        return result

    def calculate_phenotype_correlation(
        df_input: pd.DataFrame, 
        phenotype: list[float], 
        gene_sets: dict[str, set],
        threshold: float=0.0,
        is_ctrl: bool=True) -> pd.DataFrame:
        """Returns a new DataFrame needed for running GSEA and a DataFrame cleaned up of NaN corrs

        Returns:
            - df_corrs: Has the following columns
                - "accession_number" : Needed to identify the gene sets the gene belongs to
                - "corrs" : Correlation coefficient from Pearson R
                - "corrs_stats" : The p-value associated with the particular correlation coefficient
        
        """

        df_expressions = df_input.copy()
        if is_ctrl:
            expressions = df_expressions.iloc[:,4:10].to_numpy()
        else:
            expressions = df_expressions.iloc[:,10:15].to_numpy()
    
        corrs = np.zeros((np.size(expressions, 0), 2))
        for i, expression in enumerate(expressions):
            pearson_result = scipy.stats.pearsonr(phenotype, expression)
            corrs[i,:] = [pearson_result.statistic, pearson_result.pvalue]

        # Filter out NaN
        df_corrs = pd.DataFrame({
            "accession_number": df_input["accession_number"],
            "corrs": corrs[:,0],
            "corrs_stats": corrs[:,1]
        })
        is_nan = np.isnan(corrs[:,0])
        df_corrs = df_corrs[~is_nan]
        removed_genes = df_expressions.loc[is_nan, "accession_number"].values

        output_gene_sets = dict()
        for gene_set_name, gene_set in gene_sets.items():
            output_gene_sets[gene_set_name] = gene_set - set(removed_genes)

        df_expressions = df_expressions[~is_nan]

        # Applying correlation threshold filter
        # DEPRACATED
        is_low = abs(df_corrs.loc[:,"corrs"].values) < threshold
        df_corrs = df_corrs[~is_low]
        df_expressions = df_expressions[~is_low]

        return df_corrs, df_expressions, output_gene_sets

    def sort_expressions_by_rank(df_expressions, df_corrs, phenotype):
        """Sorts df_expressions based on df_corrs

        """
        # expressions = df_expressions.iloc[:,4:10]
        corrs = df_corrs.loc[:,"corrs"]
        idx_sorting = corrs.argsort()[::-1].values
        sorted_df_expressions = df_expressions.iloc[idx_sorting, :]
        sorted_df_corrs = df_corrs.iloc[idx_sorting, :]
        return sorted_df_expressions, sorted_df_corrs

    def normalize_expressions(expression: np.ndarray):
        """Normalize expression for heat map visualization

        Technically, I'm assuming that rescaling and recentering doesn't affect the correlation
        coefficient. I'm mainly doing this so that genes on the heatmap can be represented on
        a similar scale to be compared for their relative change with respect to the max deviation
        from the original expression level on Day 0.
        """

        normalized_expressions = expression.T - expression[:,0].T
        normalized_expressions = normalized_expressions / np.max(abs(normalized_expressions), axis=0)
        return normalized_expressions.T
    
    def calculate_enrichment_score(df: pd.DataFrame, gene_sets: dict[str, set], p: float) -> pd.DataFrame:
        """Calculates enrichment score for all the gene_sets

        Based on the Subramanian et al. 2005

        Arguments:
            - df_all: the dataframe that contains all the genes under study with columns
                - corr: correlation coefficients to the phenotype
                - accession_number: to check 
            - gene_sets: the name of the gene set mapped to a Python set of genes that belong to the set
            - p: exponent to control the weight of each step
        Returns:
            - df_es: Each column contains the running cumulatives enrichment score
            - max_es: numpy array of the max of each gene set
        """
        df_output = pd.DataFrame()

        for gene_set_name, gene_set in gene_sets.items():
            es_cum: float = 0
            es = [es_cum]
            nr: float = np.sum([abs(row["corrs"])**p for _, row in df.iterrows() if row["accession_number"] in gene_set])
            p_miss = 1 / (df.shape[0] - len(gene_set))
            for i, row in df.iterrows():
                if row["accession_number"] in gene_set:
                    es_cum += abs(row["corrs"])**p / nr
                else:
                    es_cum -= p_miss
                es.append(es_cum)
            df_output[gene_set_name] = es

        max_es = np.array([[gene_set_name, max(min(df_output.loc[:,gene_set_name]), max(df_output.loc[:,gene_set_name]), key=abs)] for gene_set_name, gene_set in gene_sets.items()])
        return df_output, max_es
        
    
    def calculate_rand_es_distribution(df, gene_sets, p, max_es_og, num_rand_samples: int=1000, plot_enabled: bool=False):
        # enrichment_scores = np.zeros((len(gene_sets), num_rand_samples))
        num_gene_sets = len(gene_sets)
        gene_set_names = np.empty((num_gene_sets*num_rand_samples), dtype=np.dtype("U32"))
        enrichment_scores = np.empty((num_gene_sets*num_rand_samples), dtype=np.float64)
        is_es_og_pos = [float(es[1]) > 0 for es in max_es_og]
        for i in range(num_rand_samples):
            print(f"Random ES Distribution Progress: {i} of {num_rand_samples} ({round(i/num_rand_samples*100, 2)}%)", end="\r")
            _, max_es = calculate_enrichment_score(df.sample(frac=1).reset_index(drop=True), gene_sets=gene_sets, p=p)
            is_es_pos = [float(es[1]) > 0 for es in max_es]
            for ix, (is_og_pos, is_rand_pos) in enumerate(zip(is_es_og_pos, is_es_pos)):
                if (is_og_pos != is_rand_pos):
                    max_es[ix,1] = str(-float(max_es[ix,1]))
            enrichment_scores[i*num_gene_sets:(i+1)*num_gene_sets] = max_es[:,1]
            gene_set_names[i*num_gene_sets:(i+1)*num_gene_sets] = max_es[:,0]
            sys.stdout.write("\033[K")
        
        df_rand_es_distribution = pd.DataFrame({
            "gene_set_name": gene_set_names,
            "es": enrichment_scores
        })

        # df_rand_es_distribution.to_csv("./data/processed/random_enrichment_score_distribution.csv")
        if plot_enabled:
            sns.histplot(data=df_rand_es_distribution, x="es", hue="gene_set_name", multiple="stack")
        # plt.show()
        return df_rand_es_distribution
    
    def calculate_es_p_values(max_es, df_rand_es_distribution: pd.DataFrame, gene_sets, num_genes):
        # df_pivot = df_rand_es_distribution.pivot(columns="gene_set_name", values="es")
        # print(df_pivot)
        plt.figure()
        p_values = []
        gene_set_names = []
        gene_set_ratios = []
        for i, (gene_set_name, gene_set) in enumerate(gene_sets.items()):
            df_es_gene_set = df_rand_es_distribution[df_rand_es_distribution["gene_set_name"] == gene_set_name]
            cum_count = 0
            es_values = np.sort(df_es_gene_set["es"].values)
            # print(es_values)
            counts, bin_edges = np.histogram(df_es_gene_set["es"].values, density=True)
            dx = bin_edges[1] - bin_edges[0]
            f1 = np.cumsum(counts)*dx
            es_max = float(max_es[i,1])
            is_es_og_pos = float(max_es[i,1]) > 0


            p_value_index = np.argmax(bin_edges[1:] > es_max)
            p_values.append(1 - f1[p_value_index] if is_es_og_pos else f1[p_value_index])
            gene_set_names.append(gene_set_name)
            gene_set_ratios.append(len(gene_set)/num_genes)
            # norm_cdf = scipy.stats.norm.cdf(es_values)
            # mean_es = np.mean(es_values)
            # es_max = float(max_es[i,1])
            # if es_values[-1] < es_max:
            #     p_value_index = len(es_values) - 1
            #     p_value = 1 - norm_cdf[p_value_index]
            # elif es_values[0] > es_max:
            #     p_value_index = 0
            #     p_value = norm_cdf[p_value_index]
            # elif mean_es <= es_max:
            #     p_value_index = np.argmax(es_values>float(max_es[i,1]))
            #     p_value = 1 - norm_cdf[p_value_index]
            # else:
            #     p_value_index = np.argmax(es_values<float(max_es[i,1]))
            #     p_value = norm_cdf[p_value_index]
            # # print(es_values)
            # print(p_value_index, norm_cdf[p_value_index], max_es[i,1], p_value)


            # sns.lineplot(x=es_values, y=norm_cdf)
            # plt.show()
            
            # for ix, bin_edge in enumerate(bin_edges):
            #     print(max_es[i,1])
            #     if float(bin_edge) < float(max_es[i,1]):
            #         cum_count += counts[ix]
            #     else:
            #         break
            # output.append(cum_count)

        output_df = pd.DataFrame({
            "gene_set_names": max_es[:,0],
            "es": max_es[:,1],
            "p_values": p_values,
            "gene_ratio": gene_set_ratios
        })
        print(output_df)
        # pprint.pprint(p_values)
        output_df.to_csv("./data/processed/gsea_signaling_results.csv")
        return output_df
    

    df_expressions_gene_sets = retrieve_expression_data(gene_set_paths=gene_set_paths)
    gene_sets = get_gene_sets(gene_set_paths)
    df_corrs, df_expressions, gene_sets = calculate_phenotype_correlation(df_expressions_gene_sets, 
                                                                          phenotype, gene_sets, 
                                                                          threshold=0,
                                                                          is_ctrl=is_ctrl)

    sorted_df_expressions, sorted_df_corrs = sort_expressions_by_rank(df_expressions, df_corrs, phenotype)

    if is_ctrl:
        expressions = sorted_df_expressions.iloc[:,4:10].to_numpy()
    else:
        expressions = sorted_df_expressions.iloc[:,10:15].to_numpy()
        
    ranked_list = normalize_expressions(expressions)

    df_es, max_es = calculate_enrichment_score(sorted_df_corrs, gene_sets=gene_sets, p=1)
    # df_es_random = calculate_enrichment_score(sorted_df_corrs.sample(frac=1).reset_index(drop=True), gene_sets=gene_sets, p=1)
    # plt.plot(df_es_random.loc[:,gene_set_paths.keys()])

    df_rand_es_distribution = calculate_rand_es_distribution(df_corrs, gene_sets, max_es_og=max_es, p=1)
    calculate_es_p_values(max_es, df_rand_es_distribution, gene_sets=gene_sets, num_genes=np.size(expressions,0))

    if plot_enabled:
        plt.figure("Heat map")
        sns.heatmap(ranked_list[:,1:], cmap="vlag")
        plt.title(title_heatmap)

        plt.figure("GSEA")
        plt.title(title_gsea)
        plt.plot(df_es.loc[:,gene_set_paths.keys()])
        plt.legend(gene_set_paths.keys())
        plt.show()

def dotplot():
    df_gsea = pd.read_csv("./data/processed/gsea_signaling_results.csv")
    df_gsea = df_gsea.sort_values(by="es", ascending=False)
    df_gsea["log_p_values"] = np.log(df_gsea["p_values"].values)
    df_gsea.fillna(value=-7, inplace=True)
    # sns.set_style("dark_grid")
    # sns.set_style("whitegrid")
    fig, ax = plt.subplots()
    sns.scatterplot(data=df_gsea, 
                         x="es", 
                         y="gene_set_names", 
                         size="gene_ratio", 
                         hue="log_p_values", 
                         legend="brief", 
                         sizes=(50, 300),)
    
    # ax = plt.gca()
    # ax.grid(b=True, which="major")
    plt.show()
    # fig = gsea_scatter.get_figure()
    fig.savefig("./figures/gsea_signaling_results.png", dpi=600, format="png")
    display(df_gsea)

if __name__ == "__main__":
    # compute_corr()
    # graph_corr([], "./data/processed/03_gene_sets/corr_matrix_ctrl.npy", title="Ethylene-Auxin Crosstalk Network After Ethylene-Receptor Inhibition")
    phenotype_golden_delicious_ctrl = [3, 0, 2.0349, 43.3721, 56.09, 52.0930]
    gene_set_paths = {
            "auxin_signaling": "P:auxin-activated signaling pathway", 
            "ethylene_signaling": "P:ethylene-activated signaling pathway",
            "cytokinin_signaling": "P:cytokinin-activated signaling pathway",
            "jasmonic_acid_mediation": "P:jasmonic acid mediated signaling pathway",
            "golgi_membrane_txport": "P:Golgi to plasma membrane transport"}
    is_ctrl = True
    gsea(
        phenotype=phenotype_golden_delicious_ctrl,
        gene_set_paths=gene_set_paths,
        is_ctrl=is_ctrl,
        title_heatmap="Granny Smith Control Ethylene vs Auxin",
        title_gsea="GSEA",
        plot_enabled=True)