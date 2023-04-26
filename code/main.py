from preprocess import annotate, utils, gene_sets
import logging
import datetime
if __name__ == "__main__":
    logging.basicConfig(filename="./data/test/logging.txt",
                        level=logging.ERROR)

    logging.error(f"Starting session at {datetime.datetime.now()}")
    #utils.clean_intermediates()
    input_file = "./data/processed/01_processed_fpkm_.xlsx"
    path_annotated_fpkm = "./data/processed/02_annotated_fpkm_v3.xlsx"
    properties=[
        "keywords", "go_gene_set", "protein_existence", 
        "primary_accession", "uniProtkbId", "protein", "pipeline"]
    # annotate.collect_fpkm_annotations(input_file, properties=properties, autosave_interval_count=25)
    # annotate.preprocess_fpkm_annotations(input_file, output_xlsx=path_annotated_fpkm, properties=properties)

    gene_ontologies = [
        "P:ethylene-activated signaling pathway",
        "P:auxin-activated signaling pathway"]
        
    path_go_output_dir = "./data/processed/03_gene_sets"
    gene_sets.gene_setter(input_xlsx=path_annotated_fpkm, output_dir=path_go_output_dir, gene_ontologies=gene_ontologies)