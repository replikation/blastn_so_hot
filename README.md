![logo](logo/logo.jpg)

![](https://img.shields.io/badge/nextflow-20.01.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)
![](https://img.shields.io/badge/licence-GPL--3.0-lightgrey.svg)

[![Twitter Follow](https://img.shields.io/twitter/follow/gcloudChris.svg?style=social)](https://twitter.com/gcloudChris) 

* in development use at your own risk

* this repo is a bit for the quick daily "blasts" to check up on things
  * you have a sequence(s) and reference(s) and want to check similarities?
  * either check against NCBI (might break depending on your query) or use a set of own files (recommended)
* you get a quick graphical plot in the end (see below)

```
# local
nextflow run replikation/blastn_so_hot --fasta query.fa --references references.fasta --cores 8 -profile local,docker
# NCBI query
nextflow run replikation/blastn_so_hot --fasta query.fa --cores 8 -profile local,docker
```

# example results
* blastn output gets ploted, and looks like this:
  * each "hit region" (hsp) is shown
  * chart is **not** sorted by score or best hit
![chart](logo/plot.png)


# issues
samples with the exact same name will overlay on the final plot