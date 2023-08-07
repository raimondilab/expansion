"""
Retrieve and prepare the data for the Webserver database
"""

import os
import pickle
import random
import traceback
from sqlite3 import Error
import pandas as pd
import pymysql
from gprofiler import GProfiler



# Getting data from result file
def read_data(filename):
    """
    Read data from a CSV file and return it as a Pandas DataFrame.

    Parameters:
        filename (str): The name of the CSV file to read.

    Returns:
        pd.DataFrame: Pandas DataFrame containing the data from the CSV file.
    """
    df = pd.read_csv(filename, header=0, low_memory=False, sep=",")
    return df


# Create hexadecimal colors annotation list
def create_hexadecimal_colors(protein_domains):
    """
    Create a mapping of protein domains to hexadecimal colors. Each unique protein domain
    is associated with a random color.

    Parameters:
        protein_domains (list): A list of protein domains (strings).

    Returns:
        dict: A dictionary with protein domains as keys and hexadecimal colors as values.
    """
    domain_colors = {}

    for i in range(len(protein_domains)):
        color = "%06x" % random.randint(0, 0xFFFFFF)
        protein_domain = protein_domains[i]

        if color not in domain_colors and protein_domain not in domain_colors.keys():
            domain_colors[protein_domain] = color

    return domain_colors


# Inserting Pandas DataFrames Into Databases Using INSERT
def dataframe_to_mysql(data_file):
    """
    Insert records from a Pandas DataFrame into a MySQL database table.

    Parameters:
        data_file (Pandas DataFrame): The DataFrame containing the data to be inserted.

    Returns:
        None
    """
    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user='',
                                 password='!',
                                 db='')
    # create cursor
    cursor = connection.cursor()

    # creating column list for insertion
    cols = "`,`".join([str(i) for i in data_file.columns.tolist()])

    # Insert DataFrame records one by one.
    for i, row in data_file.iterrows():
        sql = "INSERT INTO `DE` (`" + cols + "`) VALUES (" + "%s," * (len(row) - 1) + "%s)"
        cursor.execute(sql, tuple(row))
        connection.commit()


# Extract PTMs from PhosphoSite database
def extract_ptms(uniprot_accession):
    """
    Extract post-translational modifications (PTMs) from the PhosphoSite database for a given UniProt accession.

    Parameters:
        uniprot_accession (str): The UniProt accession of a protein.

    Returns:
        list: A list of PTMs for the given protein.
    """
    ptm = []

    try:
        db = "PSP.pickle"

        if os.path.exists(db):

            with open(db, 'rb') as handle:
                # Data load
                ptms = pickle.load(handle)

                if uniprot_accession in ptms:
                    ptm = ptms[str(uniprot_accession).strip()]
                    ptm = [str(k) + ":" + str(v) for k, v in ptm.items()]

    except pickle.PickleError as er:
        print(traceback.print_exc(), er)

    return ptm


# Enrichment analysis of functional (GO and other) terms
def get_enrichment_analysis(query_input):
    """
    Perform enrichment analysis of functional (Gene Ontology and other) terms using g:Profiler API.

    Parameters:
        query_input (str): The query input for which enrichment analysis will be performed.

    Returns:
        pd.DataFrame: A Pandas DataFrame with the enrichment analysis results.
    """
    enrichment = []

    try:
        gp = GProfiler(return_dataframe=True)
        enrichment = gp.profile(organism='hsapiens', no_evidences=False, query=query_input)
    except Error as er:
        print(traceback.print_exc(), er)

    return enrichment


# Pickling MySQL table Data
def sql_table_to_pickle():
    """
    Read data from a MySQL table and pickle it into a binary file.

    Parameters:
        None

    Returns:
        None
    """
    # Connect to the database
    connection = pymysql.connect(host='localhost',
                                 user='splice',
                                 password='Spliceisoforms22#@!',
                                 db='spliceform_analysis')

    df = pd.read_sql("SELECT * FROM BINDING_REGION", connection)

    # Pickling MySQL table Data
    df.to_pickle("binding_region.pkl", protocol=pickle.HIGHEST_PROTOCOL)


# Map ENSEMBL ACC TO UNIPROT ACC
def map_acc(ensemble_accession):
    """
    Map Ensembl accession to UniProt accession.

    Parameters:
        ensemble_accession (str): The Ensembl accession to be mapped.

    Returns:
        str: The corresponding UniProt accession.
    """
    uniprot_accession = None

    try:
        db = "Proteomes_ENST2UniAC.pickle"

        if os.path.exists(db):

            with open(db, 'rb') as handle:
                # Data load
                id_mapping = pickle.load(handle)

                if ensemble_accession in id_mapping:
                    uniprot_accession = id_mapping[str(ensemble_accession).strip()]

    except pickle.PickleError as er:
        print(traceback.print_exc(), er)

    return str(uniprot_accession)


# Prepare PTMs
def prepare_ptms(ptms):
    """
    Process PTMs retrieved from the PhosphoSite database and return a formatted list of PTMs
    with relevant information.

    Parameters:
        ptms (str): PTMs as a comma-separated string.

    Returns:
        list: A list of dictionaries, each containing the position, color, and type of the PTM.
    """
    processed_ptms = []

    if len(ptms) > 0:
        for ptm in ptms.replace("[", "").replace("]", "").replace("'", "").replace("'", "").split(","):
            if ptm != 'nan' and ptms[0] != "-" and len(ptm.split(":")) == 2:
                get_pos = int(ptm.split(":")[0])
                ptm_typ = ptm.split(":")[1]

                if get_pos:
                    if ptm_typ == "p":
                        processed_ptms.append({"position": get_pos, "color": "4393c3", "type": "Phosphorylation"})
                    elif ptm_typ == "ub":
                        processed_ptms.append({"position": get_pos, "color": "bf812d", "type": "Ubiquitylation"})
                    else:
                        processed_ptms.append({"position": get_pos, "color": "666666", "type": "Other"})

    return processed_ptms
