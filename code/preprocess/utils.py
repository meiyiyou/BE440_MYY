import pprint
import numpy as np
import os
import time
import logging

import Bio
from Bio import SeqIO

def clean_intermediates(
    properties=[
        "keywords", "go_gene_set", "protein_existence", "primary_accession", 
        "uniProtkbId", "protein", "pipeline", "progress"],
    intermediates_dir=f"./data/processed/02_intermediates"):

    for prop in properties:
        if os.path.exists(f"{intermediates_dir}/02_{prop}.npy"):
            os.remove(f"{intermediates_dir}/02_{prop}.npy")

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

def use_timer(func):
    def timer(*args, **kwargs):
        t_0 = time.time.now()
        func(*args, **kwargs)
        t_1 = time.time.now()
        print(f"Execution of {func.__name__!r} took {(t_1 - t_0):.4f} seconds.")
    
    return timer

if __name__ == "__main__":
    load_and_print_npy("./data/processed/02_progress.npy")