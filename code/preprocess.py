import pandas as pd
import numpy as np
import requests
import pprint
import re

def preprocess_gene_id_column(input_xlsx, output_xlsx="data/test_output.xlsx") -> None:
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

def old_code():
    gene_ref.append(m.group(1))
    m = re.search("\[gene\=(.+?)\]", line)
    gene.append(m.group(1))
    m = re.search("\[db_xref\=GeneID\:(.+?)\]", line)
    if m:
        db_xref.append(m.group(1))
    else:
        db_xref.append("None")
    m = re.search("\[protein\=(.+?)\]", line)
    protein.append(m.group(1))
    m = re.search("\[protein_id\=(.+?)\]", line)
    protein_id.append(m.group(1))
    counter += 1

    df_output = pd.DataFrame({
        "gene_reference": gene_ref,
        "gene": gene,
        "db_xref": db_xref,
        "protein": protein,
        "protein_id": protein_id
    })
    print(df_output)

def query_gene_by_xm(xm_accession: str, format: str="json") -> dict:

    r = requests.get(f"https://rest.uniprot.org/uniprotkb/search?query={xm_accession}&format={format}")
    return r["results"][0] # always returns top entry

def preprocess_gene_query(entry: dict):
    organism = entry["organism"]
    is_apple = (organism["commonName"] == "Apple") or (organism["scientificName"] == "Malus domestica")
    if is_apple:
        primaryAccession = r["results"][RESULT_NUM]["primaryAccession"]
        uniProtkbId = r["results"][RESULT_NUM]["uniProtkbId"]
        protein = r["results"][RESULT_NUM]["proteinDescription"]["recommendedName"]["fullName"]["value"]
        print(protein)

if __name__ == "__main__":
    xm1 = "XM_008378277"
    preprocess_gene_id_column("data/raw/GSE182822_Matrix_FPKM.xlsx", "data/processed/01_processed_fpkm.xlsx")
