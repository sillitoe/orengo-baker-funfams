
## Steps

Generating the list of reasonable quality FunFam alignments:

1. Select 'interesting' FunFams from the database (1 EC term, etc)

`01.test.sql`

2. Calculate summary stats (DOPS, gap percentage):

```bash
$ tail -n +2 dataset/01.test/01.test.dsv \
    | awk '{print $1, $2}' \
    | xargs -I XXX bash -c '~/git/cathpy/scripts/cath-align-summary -i /cath/data/v4_2_0/funfam/families/$(echo XXX | awk "{print \$1}")/$(echo XXX | awk "{print \$2}").reduced.sto' > dataset/01.test/01.test.alignsum.log &
```

3. Filter by summary stats:

```bash
$ tail -n +2 dataset/01.test/01.test.alignsum.log \
    | awk '{if ($4 > 80 && $5 < 50) print}' \
    | sort -k5,5nr \
    | grep -Fv '1.25.' \
    | head -n 5 \
    | awk '{print $1}' \
    | xargs -I XXX cp XXX dataset/01.test/
```

* DOPS > 80
* gap percentage < 50
* sort by decreasing gap percentage
* remove some "bin" architectures
* take top 5
* copy data file to repo

4. Create tarball

```bash
tar zcf 01.test.sto.tgz *.sto
```
