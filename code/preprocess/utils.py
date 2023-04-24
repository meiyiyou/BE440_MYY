import pprint
import numpy as np

def clean_intermediates(properties=["primaryAccession", "uniProtkbId", "protein", "proteinExistence", "progress"]):
    for prop in properties:
        if os.path.exists(f"./data/processed/02_{prop}.npy"):
            os.remove(f"./data/processed/02_{prop}.npy")

def load_and_print_npy(npy_path):
    pprint.pprint(np.load(npy_path))
    np.savetxt("./data/test/junk.csv", np.load(npy_path), delimiter=",", fmt="%s")

def parse_gb() -> list[dict]:
    gb_path = "./data/test/test.gb"
    test = Bio.GenBank.Record.Record()
    gb_record = SeqIO.read(open(gb_path, "r"), "genbank")
    print(type(gb_record.annotations), )
    pprint.pprint(gb_record.annotations)
    return [{
        "organism": gb_record.annotations["organism"],

    }]