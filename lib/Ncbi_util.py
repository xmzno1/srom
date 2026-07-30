import sys
from optparse import OptionParser
from Bio import Entrez
from Bio import Medline
from xml.etree import ElementTree
import time
# ------------------------------------
# constants
# ------------------------------------
# *Always* tell NCBI who you are
Entrez.email = "xmzno1@stu.xjtu.edu.cn"

# ------------------------------------
# Misc functions
# ------------------------------------
def search_genes(id_list,search_field=''):
    """Use ESearch to convert RefSeq or Gene symbols to standard
    Entrez IDs.
    A request to esearch.cgi is like:
    http://www.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=ID_LIST[SEARCH_FIELD]
    Return a list of Entrez IDs.
    """
    #term = " OR ".join(map(lambda x:x+"["+search_field+"]",id_list))
    term = '('+ ' OR '.join(map(lambda x:x+'[symbol]',id_list)) +') AND "human"[Organism] AND alive[prop]'
    esearch_result = Entrez.esearch(db="gene",term=term,retmod="xml")
    parsed_result = Entrez.read(esearch_result)
    return parsed_result['IdList']

def fetch_genes(id_list):
    """Fetch Entrez Gene records using Bio.Entrez, in particular epost
    (to submit the data to NCBI) and efetch to retrieve the
    information, then use Entrez.read to parse the data.
    Returns a list of parsed gene records.
    """
 
    request = Entrez.epost("gene",id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        sys.exit(-1)
 
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    efetch_result = Entrez.efetch(db="gene", webenv=webEnv, query_key = queryKey, retmode="xml")
    genes = Entrez.read(efetch_result)
    #print "Retrieved %d records for %d genes" % (len(genes),len(id_list))
    return genes

def parse_genes(genes):
    """Parse various gene information including:
    1. Species name (taxonomy name)
    2. Entrez gene ID
    3. Official symbol
    4. RefSeq IDs
    5. Offical full name
    Basically, just to go through the parsed xml data.... A big headache to figure it out...
    Return a list of dictionary.
    """
    gene_info_list = []
    for gene_data in genes:
        gene_info = {}
        # get entrez ID
        try:
            gene_info["entrez_id"] = gene_data["Entrezgene_track-info"]["Gene-track"]["Gene-track_geneid"]
        except KeyError:
            gene_info["entrez_id"] = ""
            continue
        gene_info["refseq_ids"] = []
        for comment in gene_data.get("Entrezgene_comments",[]):
            # look for refSeq annotation
            if comment.get("Gene-commentary_heading",None) == "NCBI Reference Sequences (RefSeq)":
                # get sub-comments
                for subcomment in comment.get("Gene-commentary_comment",[]):
                    for product in subcomment.get("Gene-commentary_products",[]):
                        if product.get("Gene-commentary_heading",None) == "mRNA Sequence":
                            gene_info["refseq_ids"].append(product.get("Gene-commentary_accession",""))
        # get properties
        gene_info["official_symbol"] = "" # optional
        gene_info["official_full_name"] = "" # optional
        for gene_property in gene_data.get("Entrezgene_properties",[]):
            if gene_property.get("Gene-commentary_label",None) == "Nomenclature":
                for sub_property in gene_property["Gene-commentary_properties"]:
                    if sub_property.get("Gene-commentary_label",None)  == "Official Symbol":
                        gene_info["official_symbol"] = sub_property.get("Gene-commentary_text","")
                    if sub_property.get("Gene-commentary_label",None)  == "Official Full Name":
                        gene_info["official_full_name"] = sub_property.get("Gene-commentary_text","")

        # get taxname
        try:
            gene_info["taxname"] = gene_data["Entrezgene_source"]["BioSource"]["BioSource_org"]["Org-ref"]["Org-ref_taxname"]
        except KeyError:
            gene_info["taxname"] = ""
            continue
        gene_info_list.append(gene_info)

    return gene_info_list

def print_genes (gene_info_list):
    """Print out parsed entrez gene information in tab-delimited way.
    """
    # header
    print("%s\t%s\t%s\t%s\t%s" % ("TaxonomyName","EntrezID","OfficialSymbol","RefSeqIDs","OfficialFullName"))
    for g in gene_info_list:
        print("%s\t%s\t%s\t%s\t%s" % (g["taxname"],g["entrez_id"],g["official_symbol"],",".join(g["refseq_ids"]),g["official_full_name"]))

def get_gene_symbol(entrez_id, email="xmzno1@stu.xjtu.edu.cn"):
    """
    Fetch the official gene symbol from NCBI Gene database given an Entrez Gene ID.
    
    Parameters:
        entrez_id (str or int): NCBI Entrez Gene ID
        email (str): Your email address (required by NCBI)
    
    Returns:
        str: Official gene symbol, or None if not found
    """
    # NCBI requires you to specify your email
    Entrez.email = email
    Entrez.tool = "GeneSymbolFetcher"

    try:
        # Fetch gene summary in XML format
        handle = Entrez.efetch(db="gene", id=str(entrez_id), rettype="xml")
        records = Entrez.read(handle)
        handle.close()

        # Navigate XML structure to find the gene symbol
        if records and "Entrezgene_gene" in records[0]:
            gene_ref = records[0]["Entrezgene_gene"]["Gene-ref"]
            if "Gene-ref_locus" in gene_ref:
                return gene_ref["Gene-ref_locus"]

        return None

    except Exception as e:
        print(f"Error fetching gene symbol: {e}")
        return None


# ------------------------------------
# Search pubmed
# ------------------------------------
def search_pubmed(query, max_results=10000):
    """
    Search PubMed for a given query and return a list of PMIDs.
    """
    try:
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=max_results,
            sort="relevance"
        )
        results = Entrez.read(handle)
        handle.close()
        return results.get("IdList", [])
    except Exception as e:
        print(f"Error during PubMed search: {e}")
        return []

def fetch_pubmed_details(id_list):
    """
    Fetch PubMed article details (title, authors, journal, year, abstract).
    """
    if not id_list:
        return []

    try:
        handle = Entrez.efetch(
            db="pubmed",
            id=",".join(id_list),
            rettype="medline",
            retmode="text"
        )
        records = Medline.parse(handle)
        records = list(records)
        handle.close()
        return records
    except Exception as e:
        print(f"Error fetching PubMed details: {e}")
        return []
'''---Usages---
query = "machine learning cancer"
pmids = search_pubmed(query, max_results=5)
print(f"Found {len(pmids)} articles: {pmids}")

# Respect NCBI rate limits
time.sleep(0.4)

articles = fetch_pubmed_details(pmids)
for idx, article in enumerate(articles, start=1):
    title = article.get("TI", "No title available")
    authors = ", ".join(article.get("AU", []))
    journal = article.get("JT", "Unknown journal")
    year = article.get("DP", "Unknown year")
    abstract = article.get("AB", "No abstract available")
    print(f"\nArticle {idx}:")
    print(f"Title: {title}")
    print(f"Authors: {authors}")
    print(f"Journal: {journal} ({year})")
    print(f"Abstract: {abstract[:300]}{'...' if len(abstract) > 300 else ''}")
'''