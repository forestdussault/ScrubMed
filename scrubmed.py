#!/usr/bin/env python

from Bio import Entrez
from bs4 import BeautifulSoup

import os
import time
import click
import logging
import pandas as pd

EMAIL_ADDRESS = 'forest.dussault@inspection.gc.ca'


def chunk_list(a, n):
    """
    Splits a provided list into n chunks
    :param a: List
    :param n: Number of chunks requested
    :return: A list of lists (n chunks)
    """
    k, m = divmod(len(a), n)
    try:
        chunked_list = [a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]
    except:
        chunked_list = a
    return chunked_list


def get_num_records_from_query(query, email):
    Entrez.email = email
    handle = Entrez.esearch(db='pubmed', term=query, idtype='acc')
    record = Entrez.read(handle)
    num_records = int(record['Count'])
    return num_records


def retrieve_all_pmids(query, output_directory, email):
    """
    Uses esearch to retrieve a list of PMIDs matching your search critieria
    e.g. ' "0000/01/01"[PDAT] : "3000/12/31"[PDAT] antimicrobial ' will return every available PMID on PubMed
    Note: esearch maxes out at 100,000 records, so you have to continue searching by incrementing `retstart` parameter
    """
    Entrez.email = email
    id_list = []

    # Figure out range
    num_records = get_num_records_from_query(query=query, email=email)
    if num_records == 0:
        logging.info('No records detected for query "{}". Quitting.'.format(query))
        quit()
    elif num_records > 100000:
        logging.info('Detected {} total records'.format(num_records))
    record_range = round(num_records/100000)
    if record_range == 0:
        record_range = 1

    for x in range(record_range):
        try:
            handle = Entrez.esearch(db='pubmed', retmax=100000, retstart=x, term=query,
                                    idtype='acc', tool='ScrubMed')
            record = Entrez.read(handle)
            handle.close()
            id_list.extend(record['IdList'])
            logging.info('Retrieved {} IDs from PubMed'.format(len(id_list)))
        except:
            continue

    return id_list


def efetch_record(pmid, email):
    """
    Fetches an individual record by PMID
    :param pmid: string PMID
    :param email:
    :return: record that can be parsed by bs4 to extract information like title, date, etc
    """
    Entrez.email = email
    try:
        handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml', tool='ScrubMed')
    except:
        handle = None
    return handle


def efetch_record_env(webenv, querykey, email):
    """
    Takes a webenv and query_key gathered from epost. This allows fetching a large queryset instead of an individual
    record.
    :param webenv:
    :param querykey:
    :param email:
    :return: record(s) that can be parsed with bs4
    """
    Entrez.email = email
    handle = Entrez.efetch(db='pubmed', webenv=webenv, query_key=querykey, retmode='xml', tool='ScrubMed')
    records = handle.read()
    handle.close()
    return records


def epost_records(pmid_list, email):
    """
    Takes a list of PMIDs and returns the webenv and query_key for the search that can be passed to efetch to
    actually retrieve the records
    :param pmid_list:
    :param email:
    :return webenv, query_key: pass objects to efetch
    """
    Entrez.email = email
    try:
        search_results = Entrez.read(Entrez.epost("pubmed", id=",".join(pmid_list)))
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
    except:
        webenv = None
        query_key = None
    return webenv, query_key


def parse_xml_record(xml_record):
    """
    Parses a PubMed XML formatted record and returns abstract, title, pub year, and PMID
    :param xml_record:
    :return: dictionary containing the abstract text, title, pub year, and PMID.
    """
    try:
        pmid = xml_record.findAll('pmid')[0]
        pmid = [x for x in pmid][0]
    except (IndexError, AttributeError):
        pmid = None

    try:
        abstract = xml_record.findAll('abstracttext')
        abstract = ''.join([x.text for x in abstract])
        if abstract == '':
            abstract = None
    except (IndexError, AttributeError):
        abstract = None

    try:
        title = xml_record.findAll('articletitle')
        title = ''.join([x.text for x in title])
    except (IndexError, AttributeError):
        title = None

    try:
        year = xml_record.findAll('year')[0].text
    except (IndexError, AttributeError):
        year = None

    return {'pmid': pmid, 'abstract': abstract, 'title': title, 'year': year}


def download_pipeline(output_directory, query):
    pubmed_ids = retrieve_all_pmids(query=query, output_directory=output_directory, email=EMAIL_ADDRESS)
    num_records = len(pubmed_ids)

    # Cut into n chunks
    if num_records >= 10000:
        pubmed_id_chunks = chunk_list(a=pubmed_ids, n=round(num_records/1000))
        logging.info('Split master list into {} chunks'.format(len(list(pubmed_id_chunks))))
    else:
        pubmed_id_chunks = [pubmed_ids]

    logging.info('Retrieving records via Entrez...')
    # Iterate over chunked list
    for i, chunk in enumerate(pubmed_id_chunks):
        # efetch every record in the chunk
        _webenv, _querykey = epost_records(pmid_list=chunk, email=EMAIL_ADDRESS)
        data = efetch_record_env(webenv=_webenv, querykey=_querykey, email=EMAIL_ADDRESS)

        # split the data into separate articles
        _xml_records = BeautifulSoup(data, 'lxml')
        records = _xml_records.findAll('pubmedarticle')

        # Master list containing a dictionary for each entry
        pubmed_list = []

        # Parse each of the records
        for record in records:
            # Retrieve pmid, title, abstract, publication year
            record_detail_dict = parse_xml_record(record)
            pubmed_list.append(record_detail_dict)

        # Use pandas to export everything neatly into .csv files
        df_list = []
        for dictionary in pubmed_list:
            if dictionary['abstract'] is not None:
                df = pd.DataFrame.from_records([dictionary])
                df_list.append(df)
            else:
                continue
        df = pd.concat(df_list)

        # Drop master dictionary into csv
        master_dict_filepath = os.path.join(output_directory, 'master_pubmed_dict_{}.csv'.format(i+1))
        df.to_csv(master_dict_filepath, index=False)

        # Get file stats on output file
        info = os.stat(master_dict_filepath)
        logging.info('Wrote {} kB to {}'.format((info.st_size / 1000), master_dict_filepath))


@click.command()
@click.option('-o', '--output_directory',
              required=True,
              help='Output directory to dump PubMed data')
@click.option('-q', '--query',
              required=True,
              help='Your query to send to PubMed. Note that your query must be in double quotes.'
                   ' e.g. ""0000/01/01"[PDAT] : "3000/12/31"[PDAT] antimicrobial"')
def scrubmed_cli(output_directory, query):
    # Logging
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    logging.info('Starting ScrubMed')
    start_time = time.time()

    # Validation
    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)
        logging.info('Created output directory at {}'.format(output_directory))

    # Call pipeline
    logging.info('QUERY: {}'.format(query))
    download_pipeline(output_directory=output_directory, query=query)

    # Completion
    logging.info('Total time elapsed: {0:.{1}f} seconds'.format(time.time() - start_time, 2))


if __name__ == '__main__':
    scrubmed_cli()
