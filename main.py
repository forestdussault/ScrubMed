from Bio import Entrez
from bs4 import BeautifulSoup

import os
import sys
import time
import logging
import pandas as pd
import _pickle as pickle

# Prevent pickle errors resulting from pickling large objects
sys.setrecursionlimit(1000000)


def load_pickle_list(pickle_file):
    logging.info('Loading {}...'.format(pickle_file))
    f = open(pickle_file, 'rb')
    pickle_list = pickle.load(f)
    f.close()
    logging.info('Loaded pickle successfully')
    return pickle_list


def dump_to_pickle(filename, to_pickle):
    pickled_object_location = filename
    with open(pickled_object_location, 'wb') as file:
        pickle.dump(to_pickle, file)
    logging.info('Successfully dumped {} to {}'.format(to_pickle, filename))


def chunk_list(a, n):
    k, m = divmod(len(a), n)
    return [a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]


def retrieve_all_pmids(searchterm):
    """
    Uses esearch to retrieve a list of PMIDs matching your search critieria
    e.g. ' "0000/01/01"[PDAT] : "3000/12/31"[PDAT] antimicrobial ' will return every available PMID on PubMed
    """
    Entrez.email = 'forest.dussault@inspection.gc.ca'
    # ESearch for every single pubmed record: "0000/01/01"[PDAT] : "3000/12/31"[PDAT]
    # esearch maxes out at 100,000 records, so you have to continue searching by incrementing `retstart`
    id_list = []
    for x in range(16):
        try:
            handle = Entrez.esearch(db='pubmed', retmax=100000, retstart=x, term=searchterm, idtype='acc')
            record = Entrez.read(handle)
            handle.close()
            id_list.extend(record['IdList'])
            logging.info('Retrieved {} records from PubMed'.format(len(id_list)))
        except:
            continue

    # Dump to pickle to avoid doing this again
    dump_to_pickle(filename=os.path.join('/home/forest/PycharmProjects/ScrubMed', 'all_pubmed_ids.pickle'),
                   to_pickle=id_list)


def efetch_record(pmid):
    """
    Fetches an individual record by PMID
    :param pmid: string PMID
    :return: record that can be parsed by bs4 to extract information like title, date, etc
    """
    Entrez.email = 'forest.dussault@inspection.gc.ca'
    try:
        handle = Entrez.efetch(db='pubmed', id=pmid, retmode='xml', rettype='abstract')
    except:
        handle = None
    return handle


def efetch_record_env(webenv, querykey):
    """
    Takes a webenv and query_key gathered from epost. This allows fetching a large queryset instead of an individual
    record.
    :param webenv:
    :param querykey:
    :return: record(s) that can be parsed with bs4
    """
    Entrez.email = 'forest.dussault@inspection.gc.ca'
    handle = Entrez.efetch(db='pubmed', webenv=webenv, query_key=querykey, retmode='xml')
    records = handle.read()
    handle.close()
    return records


def epost_records(pmid_list):
    """
    Takes a list of PMIDs and returns the webenv and query_key for the search that can be passed to efetch to
    actually retrieve the records
    :param pmid_list:
    :return webenv, query_key: pass objects to efetch
    """
    Entrez.email = 'forest.dussault@inspection.gc.ca'
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
        abstract = xml_record.findAll('abstracttext')[0]
        abstract = [x for x in abstract][0]
    except (IndexError, AttributeError):
        abstract = None

    try:
        title = xml_record.findAll('articletitle')[0]
        title = [x for x in title][0]
    except (IndexError, AttributeError):
        title = None

    try:
        year = xml_record.articledate.findAll('year')[0]
        year = [x for x in year][0]
        if year == '':
            year = xml_record.pubmedpubdate.year.text
    except (IndexError, AttributeError):
        year = None

    return {'pmid': pmid, 'abstract': abstract, 'title': title, 'year': year}


def main():
    # retrieve_all_pmids('"0000/01/01"[PDAT] : "3000/12/31"[PDAT] antimicrobial')

    pubmed_ids = load_pickle_list('all_pubmed_ids.pickle')
    logging.info('Total records available: {}'.format(len(pubmed_ids)))

    # Cut into n chunks
    pubmed_id_chunks = chunk_list(a=pubmed_ids, n=1000)
    logging.info('Split master list into {} chunks'.format(len(list(pubmed_id_chunks))))

    # Iterate over chunked list
    for i, chunk in enumerate(pubmed_id_chunks):

        logging.info('{}. Retrieving {} records via Entrez...'.format(i, len(chunk)))

        # efetch every record in the chunk
        _webenv, _querykey = epost_records(chunk)
        data = efetch_record_env(_webenv, _querykey)

        # split the data into separate articles
        _xml_records = BeautifulSoup(data, 'lxml')
        records = _xml_records.findAll('pubmedarticle')

        # Master list containing a dictionary for each entry
        pubmed_list = []
        # parse each of the records
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
        master_dict_filepath = os.path.join('/mnt/nas2/users/forest/pubmed_abstract_scrape',
                                            'master_pubmed_dict_{}.csv'.format(i))
        df.to_csv(master_dict_filepath, index=False)

        # Get file stats on pickle
        info = os.stat(master_dict_filepath)
        logging.info('Filesize: {} kB'.format(info.st_size/1000))


if __name__ == '__main__':
    logging.basicConfig(format='\033[92m \033[1m %(asctime)s \033[0m %(message)s ',
                        level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
    start_time = time.time()

    main()

    logging.info('Total time elapsed: {0:.{1}f} seconds'.format(time.time() - start_time, 2))
