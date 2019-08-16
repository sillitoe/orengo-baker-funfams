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
