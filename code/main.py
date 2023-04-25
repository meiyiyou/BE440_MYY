from preprocess import annotate, utils
import logging
import datetime
if __name__ == "__main__":
    logging.basicConfig(filename="./data/test/logging.txt",
                        level=logging.ERROR)

    logging.error(f"Starting session at {datetime.datetime.now()}")
    #utils.clean_intermediates()
    input_file = "./data/processed/01_processed_fpkm.xlsx"
    output_file = "./data/processed/02_annotated_fpkm.xlsx"
    properties=["keywords", "go_gene_set", "protein_existence", "primary_accession", "uniProtkbId", "protein", "pipeline"]
    annotate.collect_fpkm_annotations(input_file, properties=properties, autosave_interval_count=25)
    # annotate.preprocess_fpkm_annotations(input_file, output_xlsx=output_file, properties=properties)