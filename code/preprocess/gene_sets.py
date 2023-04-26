"""This module generates in-house gene sets based on ontologies.

By creating gene sets based on the the GO:0009755 (hormone-mediated signaling pathway),
we can identify if there are any entities that are part of multiple hormone-mediated 
pathways. These are our "cross-talk regions" of interest (de Anda-Jáuregui et al. 
BMC Systems Biology).

The output of this module can be used for analysis
    -   Generated a network based on gene ontology (nodes are individual genes positioned
        based on the particular ontology node)
    -   Genes that are part of multiple ontologies can be an edge between the ontology node
        Ontologies also have defined relations/edges already between themselves
    -   GSEA is possible from these gene sets
"""
import logging
import pandas as pd
import numpy as np
import re

def gene_setter(input_xlsx, output_dir, gene_ontologies) -> None:
    """
    Args:
        input_xlsx: Annotated xlsx output from annotate.py
        output_dir: WITHOUT TRAILING SLASH! RELATIVE TO PROJECT ROOT DIR.
        gene_ontologies: List of strings that correspond to gene ontologies that you want sets of

    Returns:
        None. Will save a xlsx file that only has gene of a particular set.
    """
    # STEP 1: 
    df_total = pd.read_excel(input_xlsx, sheet_name="GrannySmith_GS")
    for go in gene_ontologies:
        NUM_ENTRIES = np.size(df_total,0)
        gene_set_members = np.empty((NUM_ENTRIES,1), dtype=np.dtype("U32"))
        for i, entry in df_total.iterrows():
            go_gene_set = entry["go_gene_set"]
            if type(go_gene_set) is str:
                # if re.search(go, go_gene_set):
                #     # print(go)
                #     # print(entry["go_gene_set"])
                if re.search(go, go_gene_set):
                    # print(entry["accession_number"])
                    gene_set_members[i] = entry["accession_number"]

            if i % 250 == 1:
                print(f"{go}: {round(i / NUM_ENTRIES),2}%...", end="\r")

        output = [gene for gene in gene_set_members[:,0] if gene != ""]
        with open(f"{output_dir}/{go}.txt", "w+") as f:
            f.write("\n".join(output))

if __name__ == "__main__":
    gene_set_logger = logging.Logger(name="./data/test/gene_set_logger.txt")
    gene_set_logger.setLevel(logging.INFO)
    input_xlsx = "./data/processed/02_annotated_fpkm_v3.xlsx"
    output_dir = "./data/test/03_gene_sets"
    gene_ontologies = ["P:auxin-activated signaling pathway", "P:ethylene-activated signaling pathway"]
    gene_setter(input_xlsx, output_dir, gene_ontologies)




