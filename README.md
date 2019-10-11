# Orengo/Baker FunFam collaboration

This repo is currently here to store data for a working collaboration between the [Baker group](https://www.bakerlab.org/)
and the [Orengo group](http://orengogroup.info/) regarding modelling protein structures for Functional Family (FunFam)
alignments in CATH.

## Datasets

### 01.test

Initial set of 5 FunFam protein sequence alignments (STOCKHOLM format) to use for general testing: 

[`dataset/01.test/01.test.sto.tgz`](https://github.com/sillitoe/orengo-baker-funfams/blob/master/dataset/01.test/01.test.sto.tgz)

Selection criteria:

* 1 EC annotation
* reasonable number of sequences
* reasonable sequence diversity
* not too high gap percentage

### 02.medium

Wider set of FunFam alignments to test different aspects of pipeline:

Large alignments

Selection criteria:

* 1 EC annotation
* reasonable number of sequences
* reasonable sequence diversity
* not too high gap percentage

Large sequences

Plus:

* One example taken from each CATH fold (with above criteria)

Setup local python environment:

```sh
python3 -m venv venv
. venv/bin/activate
pip install -r requirements.txt
```

Run script to generate testset:

```sh
$ PYTHONPATH=. python3 scripts/02.medium.py

LOADING FUNFAM DATA
-------------------
Loading cached Funfam data from 'cath.v4_2_0.funfam.json' ...
Found 68065 funfams from 2623 unique superfamilies and 1327 unique folds

APPLYING GLOBAL FILTERS
-----------------------
High sequence diversity (seed_dops_score >= 70):             19053 (from 68065)

CREATING DATASET CATEGORIES
---------------------------
Small Funfams (num_members < 50):                            2633 (from 19053)
Large Funfams (num_members >= 5000):                         1492 (from 19053)
Funfams with structure (rep_source_id="cath"):               8464 (from 19053)
Funfams without structure (rep_source_id="uniprot"):         10589 (from 19053)

WRITING FUNFAM ALIGNMENTS
-------------------------

Working on 2633 funfams with low sequences ...
----------------------------------------------
Funfams with unique topologies:                              476 (from 1327)
GET https://www.cathdb.info/version/v4_2_0/api/rest/superfamily/1.10.10.10/funfam/20119/files/seed_alignment
Saving Funfam alignment 1.10.10.10-ff-20119
[--8<--]

Working on 1492 funfams with high sequences ...
-----------------------------------------------
Funfams with unique topologies:                              178 (from 1327)
GET https://www.cathdb.info/version/v4_2_0/api/rest/superfamily/1.10.1030.10/funfam/3118/files/seed_alignment
[--8<--]

Working on 8464 funfams with structure ...
------------------------------------------
Funfams with unique topologies:                              453 (from 1327)
GET https://www.cathdb.info/version/v4_2_0/api/rest/superfamily/1.10.1000.11/funfam/894/files/seed_alignment
[--8<--]

Working on 10589 funfams with no structure ...
----------------------------------------------
Funfams with unique topologies:                              42 (from 1327)
GET https://www.cathdb.info/version/v4_2_0/api/rest/superfamily/1.10.1510.10/funfam/1175/files/seed_alignment
[--8<--]

WRITING DATASETS
----------------
Total number of funfams                                      1149 (from 68065)
Writing TSV datafile 'dataset/02.medium/02.medium.all.tsv'
DONE

```

## Notes

### Alignment format

These files should adhere to the general [Stockholm format](https://en.wikipedia.org/wiki/Stockholm_format) rules, so you should be able to use any parser to convert them to FASTA.

We found that the alignment model within [BioPython](https://biopython.org/wiki/AlignIO) didn't provide the flexibility we needed so we added an alignment module to [`cathpy`](https://github.com/UCL/cathpy) that might be useful (see [`cathpy.align`](https://cathpy.readthedocs.io/en/latest/align.html)).

Converting from STOCKHOLM to FASTA with this package:

```python
from cathpy.align import Align
aln = Align.from_stockholm('path/to/aln.sto')
aln.write_fasta('path/to/aln.fa')
```
