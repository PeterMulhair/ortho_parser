# Collection of scripts for gene family setup and analyses

## Starting dataset

1. Download a set of proteomes

2. Use OrthoFinders [primary_transcript.py](https://github.com/davidemms/OrthoFinder/blob/master/tools/primary_transcript.py) script to get longest isoforms for each proteome:

```
for i in *fa; do; python primary_transcript.py $i; done
```

This outputs to a directory called `primary_transcripts/`

3. Copy the `rename_genes.py` script from this repo to the `primary_transcripts/` directory and run using `python rename_genes.py`

This outputs simplified renamed gene fasta files.

## Running OrthoFinder

Using the output of the `rename_genes.py` script (i.e. `primary_transcripts/genes_renamed/`) run OrthoFinder in whatever way you see fit.

## Building gene trees from single copy orthologs

`infer_trees.py` is a script to build alignments and trees for a given set of gene famillies.

This script requires some input:
```
python infer_trees.py --input <path to directory containing gene family fasta files> --output <path to directory you want output to go> --threads <number of threads to run in parallel>
```

It's recommended to run on multiple threads to speed up the pipeline, espcially if you have many gene families.