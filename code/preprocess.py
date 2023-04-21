import pandas as pd
import numpy as np
import requests
import pprint
import re
import logging
import os
import pathlib
import sys
import time

def preprocess_gene_id_column(input_xlsx, output_xlsx="data/test_output.xlsx") -> None:
    """Update original raw FPKM data to split the GENE_ID column
    """
    with pd.ExcelWriter(output_xlsx) as writer:
        for sheet in ["GrannySmith_GS", "GoldenDelicious_GD", "Fuji_Fj"]:
            df = pd.read_excel(input_xlsx, sheet_name=sheet)
            gi = np.empty(np.size(df,0), dtype=np.dtype("U32"))
            database = np.empty(np.size(df, 0), dtype=np.dtype("U32"))
            accession_number = np.empty(np.size(df, 0), dtype=np.dtype("U32"))
            version = np.empty(np.size(df, 0))

            for i, row in enumerate(df.loc[:, "GENE_ID"]):
                gi[i] = re.search("gi\|(\d+?)\|", row).group(1)
                database[i] = re.search("gi\|\d+\|(.+?)\|", row).group(1)
                accession_number[i] = re.search("gi\|\d+\|\w+\|(\w+?),", row).group(1)
                version[i] = re.search("gi\|\d+\|\w+\|\w+?,(.+?)\|", row).group(1)

            df_output = pd.DataFrame({
                "gi": gi,
                "database": database,
                "accession_number": accession_number,
                "accession_version": version,    
            })

            df_output = df_output.join(df.iloc[:,1:])
            df_output.to_excel(writer, sheet_name=sheet, index=False)

def query_gene_by_xm(xm_accession: str, format: str="json") -> dict | None:
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={xm_accession}&format={format}").json()
    for result in r["results"]:
        organism = result["organism"]
        # Check if entry is apple 
        # TODO: Can include other filters and checks later too
        try:
            is_apple = (organism["commonName"] == "Apple") or (organism["scientificName"] == "Malus domestica")
            if is_apple:
                return result
        except Exception as err:
            logging.error(err)
            logging.error(f"accession: {xm_accession}")

def get_data_addons(input_df):
    NUM_ENTRIES = np.size(input_df, 0)

    def get_protein(entry):
        try:
            return entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
        except KeyError as key_err:
            return entry["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
    
    def load_values(prop, num_entries=NUM_ENTRIES, data_type="U32"):
        property_path = f"./data/processed/02_{prop}.npy"
        if os.path.exists(property_path):
            return np.load(property_path)
        else:
            with open(property_path, "w") as f:
                pass
            return np.empty(num_entries, dtype=np.dtype(data_type))
    
    data_addons = {
        "primaryAccession": {
            "getter": lambda x: x["primaryAccession"],
            "values": load_values("primaryAccession")
        },
        "uniProtkbId": {
            "getter": lambda x: x["uniProtkbId"],
            "values": load_values("uniProtkbId")
        },
        "protein": {
            "getter": get_protein,
            "values": load_values("protein", data_type="U128")
        },
        "proteinExistence": {
            "getter": lambda x: x["proteinExistence"],
            "values": load_values("proteinExistence", )
        }
    }

    return data_addons

def update_data_addons(data_addons):
    for prop in data_addons.keys():
        np.save(f"./data/processed/02_{prop}.npy", data_addons[prop]["values"])

def get_progress_index(progress_file_path):
    if os.path.exists(progress_file_path):
        return np.load(progress_file_path)[0]
    else:
        np.save(progress_file_path, np.array([0]), allow_pickle=True)
        return 0

def update_progress_index(progress_index, progress_file_path):
    np.save(progress_file_path, np.array([progress_index]))

def clean_intermediates():
    properties = ["primaryAccession", "uniProtkbId", "protein", "proteinExistence", "progress"]
    for prop in properties:
        if os.path.exists(f"./data/processed/02_{prop}.npy"):
            os.remove(f"./data/processed/02_{prop}.npy")

def collect_fpkm_annotations(input_xlsx,
                            progress_file_path="./data/processed/02_progress.npy", 
                            autosave_interval_count=5) -> None:
    """Update table with annotations of proteins and genes associated with FPKM mRNA readout"""
    
    df = pd.read_excel(input_xlsx, sheet_name="GrannySmith_GS")
    NUM_MRNA = np.size(df,0)

    data_addons = get_data_addons(df)
    progress_index = get_progress_index(progress_file_path)
    for i, accession_number in enumerate(df.loc[progress_index:, "accession_number"], start=progress_index):
        print(f"Progress: Gene {accession_number} ({i} of {NUM_MRNA}) = {round(i/NUM_MRNA*100, 2)}%", end="\r")
        gene_entry = query_gene_by_xm(accession_number)
        try:
            properties = ["primaryAccession", "uniProtkbId", "protein", "proteinExistence"]
            if gene_entry is not None:
                for prop in properties:
                    data_addons[prop]["values"][i] = data_addons[prop]["getter"](gene_entry)
            else:
                for prop in properties:
                    data_addons[prop]["values"][i] = None
        except Exception as err:
            logging.error(f"{err} for gene_accession {accession_number}")
        
        if i % autosave_interval_count == 1:
            update_data_addons(data_addons)
            update_progress_index(i, progress_file_path)
        
        sys.stdout.write("\033[K") #clear line 
        
    
    update_data_addons(data_addons)
    update_progress_index(0, progress_file_path)
    print("\033[92m Done!\033[0m")

def load_and_print_npy(npy_path):
    pprint.pprint(np.load(npy_path))
    np.savetxt("./data/test/junk.csv", np.load(npy_path), delimiter=",", fmt="%s")

def update_fpkm_annotations(input_xlsx,
                            properties,
                            output_xlsx="./data/processed/test_02_annotated_fpkm.xlsx"):

    df_input = pd.read_excel(input_xlsx)
    data_addons = get_data_addons(df_input)
    df_addons = pd.DataFrame({
        prop: data_addons[prop]["values"] for prop in properties
    })

    with pd.ExcelWriter(output_xlsx) as writer:
        for sheet in ["GrannySmith_GS", "GoldenDelicious_GD", "Fuji_Fj"]:
            df_output = df_input.join(df_addons)
            df_output.to_excel(writer, sheet_name=sheet, index=False)

if __name__ == "__main__":
    logging.basicConfig(filename="./data/test/logging.txt",
                        level=logging.ERROR)
    logging.message(f"Starting preprocessing session {time.time()}")
    collect_fpkm_annotations("./data/processed/01_processed_fpkm.xlsx")