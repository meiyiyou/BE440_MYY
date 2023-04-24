"""This module is responsible fetching data from online databases.

The fetch.py module is responsible for fetching data from online databases like
NCBI (via eutils) and UniProt. This module also does some light preprocessing.

The output of the fetch module should be a FetchedData object, which is formatted
to play nicely with other modules in this pipeline, namely the annotate.py module.
The motivation for this FetchedData class is to standardize the presentation of
data retrieved from various database formats (NCBI and UniProt present information
differently, in different formats, and different amounts as well).

Although minimal query parameters need to be passed into fetch.py many keyword arguments
are provided for fine-tuned queries.

Typical usage example:

    fetch.get_fetched_data(
        accession_number, db)

"""
import requests
import urllib
from enum import Enum
from collections import UserDict
from abc import ABC, abstractmethod
import logging
import pprint
import re

import Bio.GenBank.Record
from Bio import SeqIO

# TODO: Consider turning this into a dataclass
class FetchedData(dict):
    """Data container for standardized data format for the apple project.
            "keywords": from uniprot keywords, if available
            "go_gene_set"
            "protein_existence"
            "primary_accession"
            "uniProtkbId"
            "protein"
    """

class Database(Enum):
    NCBI = 1
    UNIPROT = 2
    GENE_ONTOLOGY = 3

def get_fetched_data(accession_number, db: Database=Database.UNIPROT) -> FetchedData | None:
    """ Assigns the data request to the appropriate preprocessor
    """
    match db:
        case Database.NCBI:
            processor = NcbiPreprocessor()
        case Database.UNIPROT:
            processor = UniprotPreprocessor()
        case Database.GENE_ONTOLOGY:
            processor = GeneOntologyPreprocessor()
        
    result = processor.run(accession_number)
    if result is None:
        result = NcbiPreprocessor().run(accession_number)
    return result

class Preprocessor(ABC):
    """Abstract Base Class for the preprocessor associated with each database

    The preprocessor for each database is responsible for three key steps:

        1. Loading the data source
        2. (optional) Checking if it's an apple (looser restrictions can be implemented)
        3. Formatting the raw data into a FetchedData format
    """

    @abstractmethod
    def _01_load_data_source(
        self,
        accession_number: str, 
        url_base: str) -> dict | None:
        raise NotImplementedError
    
    # @abstractmethod
    # def _02_check_if_apple(downloaded_data):
    #     raise NotImplementedError

    @abstractmethod
    def _03_format_fetched_data(self, raw_data: dict) -> FetchedData | None:
        """Directly formats the output of _01_load_data_source
        """
        raise NotImplementedError()

    def _get_gene_sets(self, entry) -> tuple:
        """Returns (keywords, go_gene_set) tuple given Uniprot JSON entry
        entry must be a uniprot entry returned in JSON format
        """

        ontologies = {
            "Molecular function": "F", 
            "Biological process": "P", 
            "Cellular component": "C"}

        keywords = set()
        go_gene_set = set()
        for keyword in entry["keywords"]:
            if (cat := keyword["category"]) in ontologies.keys():
                keywords.add(f"{ontologies[cat]}:{keyword['name']}")

        for x_ref in entry["uniProtKBCrossReferences"]:
            if x_ref["database"] == "GO":
                for prop in x_ref["properties"]:
                    if prop["key"] == "GoTerm":
                        go_gene_set.add(prop["value"])

        return keywords, go_gene_set # TODO: Consider making this a dict that later gets unpacked
    
    def run(self, accession_number: str) -> FetchedData | None:
        raw_data = self._01_load_data_source(accession_number)
        fetched_data = self._03_format_fetched_data(raw_data) if raw_data else None
        return fetched_data

class NcbiPreprocessor(Preprocessor):
    def _01_load_data_source(
            self,
            accession_number: str, 
            url_base="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            rettype="gb"):
        
        SUB_DB = "nuccore"
        RETMODE = "text"

        gb_path = f"./data/test/{accession_number}.gb"
        r = None
        while r is None:
            try:
                r = urllib.request.urlretrieve(f"{url_base}?db={SUB_DB}&id={accession_number}&rettype={rettype}&retmode={RETMODE}", gb_path)
            except Exception as err:
                logging.error(f"{err}")

        record: Bio.SeqRecord.SeqRecord = SeqIO.read(gb_path, "genbank")
        features: list[Bio.SeqFeature.SeqFeature] = record.features

        cds_feature = None
        for feature in features:
            if feature.type == "CDS":
                cds_feature = feature
        

        return cds_feature.qualifiers if cds_feature else None

    def _03_format_fetched_data(self, raw_data):
        search_term = re.sub(r"\[|\]|\:", "", raw_data["product"][0])
        uni_data = self._get_go_terms(search_term)
        try: 
            formatted_data = FetchedData({
                "protein": raw_data["product"][0],
                "pipeline": "ncbi_to_uniprot",
                **uni_data
            })
        except KeyError:
            logging.error(f"KeyError: The dict output of _01_load_data_source did not have all the required fields")

        return formatted_data

    def _get_go_terms(self, product_search_term) -> dict | None:
        r = None
        url_base="https://rest.uniprot.org/uniprotkb/search?query="
        while r is None:
            try:
                print(f"{url_base}{product_search_term}&format=json")
                r = requests.get(f"{url_base}{product_search_term}&format=json").json()
                entries = r["results"] # Multiple results can be returned for a given query
                entry = entries[0] if len(entries) > 0 else None
            except ConnectionError as err:
                logging.error(f"{err}. product_search_term: {product_search_term}")

        if entry:
            keywords, go_gene_set = self._get_gene_sets(entry) 
            uni_data = {
                "keywords": keywords,
                "go_gene_set": go_gene_set,
                "protein_existence": entry["proteinExistence"],
                "primary_accession": entry["primaryAccession"],
                "uniProtkbId": entry["uniProtkbId"]
            }
            return uni_data
        else:
            return {}
            
        
class UniprotPreprocessor(Preprocessor):
    def _01_load_data_source(
            self,
            accession_number, 
            url_base="https://rest.uniprot.org/uniprotkb/search?query=",
            rettype="json") -> dict | None:
        r = None
        while r is None:
            try:
                r = requests.get(f"{url_base}{accession_number} AND (organism_id:3750)&format={rettype}").json()
                entries = r["results"] # Multiple results can be returned for a given query
                return entries[0] if len(entries) > 0 else None
            except ConnectionError as err:
                logging.error(f"{err}. accession: {accession_number}")

    def _get_protein(self, entry):
        try:
            return entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
        except KeyError as key_err:
            return entry["proteinDescription"]["submissionNames"][0]["fullName"]["value"]

    def _03_format_fetched_data(self, raw_data):
        keywords, go_gene_set = self._get_gene_sets(raw_data)
        formatted_data = FetchedData({
            "keywords": keywords,
            "go_gene_set": go_gene_set,
            "protein_existence": raw_data["proteinExistence"],
            "primary_accession": raw_data["primaryAccession"],
            "uniProtkbId": raw_data["uniProtkbId"],
            "protein": self._get_protein(raw_data),
            "pipeline": "uniprot",
        })

        return formatted_data

def _preprocess_ncbi() -> FetchedData:
    ### Step 1: Download

    ### Step 2: Check if apple

    ### Step 3
    fetch.download_gb(accession_number)
    results = parse_gb()
    out_fetched_data = FetchedData()
    pass

if __name__ == "__main__":
    logging.basicConfig(filename="./data/test/logging.txt",
                        level=logging.INFO)

    test_accessions = ["XM_008389256", "XM_008368189", "XM_008360699"]
    for test in test_accessions:
        pprint.pprint(get_fetched_data(test))