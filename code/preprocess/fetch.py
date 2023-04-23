def download_gb(accession_number) -> None:
    url_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    db = "nuccore"
    rettype = "gb"
    retmode = "text"
    urllib.request.urlretrieve(f"{url_base}?db={db}&id={accession_number}&rettype={rettype}&retmode={retmode}", "./data/test/test.gb")

def query_gene_by_xm(xm_accession: str, db: str="uniprot", format: str="json") -> dict | None:
    def _fetch_data(db=db):
        r = None
        while r is None:
            try:
                match db:
                    case "uniprot":
                        base_url = "https://rest.uniprot.org/uniprotkb/search?query="
                        r = requests.get(f"{base_url}{xm_accession}&format={format}").json()
                        results = r["results"] # Multiple results can be returned for a given query
                    case "ncbi":
                        fetch.download_gb(xm_accession)
                        results = parse_gb()
            except ConnectionError as err:
                logging.error(f"{err}. accession: {xm_accession}")
        return r
    
    def _check_apple(result: dict):
        match 
    for result in results:
        is_apple = _check_apple(result)
        organism = result["organism"]
        try:
            # Check if entry is apple 
            # BUG: KeyError can be a problem here
            is_apple = (organism["commonName"] == "Apple") or (organism["scientificName"] == "Malus domestica")
            if is_apple:
                return result
        except Exception as err:
            logging.error(err)
            logging.error(f"accession: {xm_accession}")