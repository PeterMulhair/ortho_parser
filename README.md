# ortho_parser
Scripts to parse OrthoFinder output

## Usage

```
usage: python ortho_parser.py [-h] -i INPUT -t TREE [-s SPECIES] [-r ROOT [ROOT ...]]
                       [-o OUTPUT] [-d THREADS]

required arguments:
  -i INPUT, --input INPUT
                        Orthofinder output directory
  -t TREE, --tree TREE  Input species tree

optional arguments:
  -h, --help            show this help message and exit
  -s SPECIES, --species SPECIES
                        Name of species of interest to infer gene stats on
  -r ROOT [ROOT ...], --root ROOT [ROOT ...]
                        Where to root unrooted species tree
  -o OUTPUT, --output OUTPUT
                        Name of output dir
  -d THREADS, --threads THREADS
                        Number of threads to run in parallel
```