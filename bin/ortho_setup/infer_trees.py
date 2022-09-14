import os
import glob
import argparse
from subprocess import call as unix
from joblib import Parallel, delayed

parse = argparse.ArgumentParser()

parse.add_argument("-i", "--input",type=str, help="SCO directory from OrthoFinder output",required=True)
parse.add_argument("-o", "--output",type=str, help="name of output directory",required=True)
parse.add_argument("-t", "--threads",type=str, help="Number of threads to run in parallel",required=True)

args = parse.parse_args()

threads = int(args.threads)

os.makedirs(args.output, exist_ok=True)
if args.output[-1] != '/':
    args.output = args.output + '/'

os.makedirs(args.output + 'aligned_genes', exist_ok=True)
os.makedirs(args.output + 'gene_trees', exist_ok=True)
    
if args.input[-1] != '/':
    args.input = args.input + '/'

fasta_list = []
for fasta in glob.glob(args.input + '*fa'):
    fasta_list.append(fasta)

def align(fasta):
    fa = fasta.split('/')[-1]
    OG = fa.split('.')[0]
    unix("mafft --quiet --auto " + fasta + " > " + args.output + "aligned_genes/" + OG + ".aln", shell=True)

def iqtree(fasta):
    OG = fasta.split('/')[-1].split('.')[0]
    unix('iqtree -s ' + args.output + 'aligned_genes/' + OG + '.aln -B 1000 --prefix ' + args.output + 'gene_trees/' + OG, shell=True)
    

Parallel(n_jobs=threads)(delayed(align)(fa) for fa in fasta_list)   
Parallel(n_jobs=threads)(delayed(iqtree)(fa) for fa in fasta_list)
