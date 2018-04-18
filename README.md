## ScrubMed
This is a simple command line utility that makes use of BioPython to
retrieve a set of abstracts, titles, PMIDs, and publication years for
a given query. The results are dumped into manageble *.csv* files.

Keep in mind that it is possible to retrieve a very large number of
records if your query is broad.

### Guidelines
Read this: https://www.ncbi.nlm.nih.gov/books/NBK25497/

**Limit large jobs to either weekends or between 9:00 PM and 5:00 AM
Eastern time during weekdays.**

### Installation
```
pip install ScrubMed
```

### Usage
```
Usage: scrubmed.py [OPTIONS]

Options:
  -o, --output_directory PATH  Output directory to dump PubMed data
                               [required]
  -q, --query TEXT             Your query to send to PubMed. Note that your
                               query must be in double quotes. e.g.
                               ""0000/01/01"[PDAT] : "3000/12/31"[PDAT]
                               antimicrobial"  [required]
  --help                       Show this message and exit.
```

#### Example usage
This will retrieve all records from 2018 with the keywords `"antimicrobial resistance"`.
All of the resulting output will be dumped into the specified output folder.
```
scrubmed.py -o /home/user/scrubmed_output -q ""2018/01/01"[PDAT] : "2018/12/31"[PDAT] antimicrobial resistance"
```