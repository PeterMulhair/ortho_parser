#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import toytree
import toyplot
import toyplot.svg
import toyplot.pdf
import numpy as np
from Bio import SeqIO
from ete3 import Tree, PhyloTree
from subprocess import call as unix
from joblib import Parallel, delayed
from collections import defaultdict, Counter
from subprocess import Popen, PIPE, check_output

# Author: Peter Mulhair
# Date: 17/08/2022
# Usage: python3

'''
Script to parse OrthoFinder output, annotate 
gene gain and loss across the species tree & 
get stats significant to a species of interest
'''

parse = argparse.ArgumentParser()

parse.add_argument("-i", "--input",type=str, help="Orthofinder output directory",required=True)
parse.add_argument("-t", "--tree",type=str, help="Input species tree",required=True)
parse.add_argument("-s", "--species",type=str, help="Name of species of interest",required=False)
parse.add_argument("-r", "--root",nargs="+", help="Where to root unrooted species tree",required=False)
parse.add_argument("-o", "--output",type=str, help="Name of output dir",default="output")
parse.add_argument("-d", "--threads",type=str, help="Number of threads to run in parallel",required=False)

args = parse.parse_args()

if args.threads:
    threads = int(args.threads)

os.makedirs(args.output, exist_ok=True)
if args.output[-1] != '/':
    args.output = args.output + '/'

if args.input[-1] != '/':
    args.input = args.input + '/'

    
#If required, root tree on species given
if args.root:
    unrooted_tree = Tree(args.tree, format = 1)
    if len(args.root) == 1:
        unrooted_tree.set_outgroup(args.root[0])
        unrooted_tree.write(format=1, outfile=args.output + "rooted_tree.nwk")
        input_tree = Tree(args.output + "rooted_tree.nwk", format = 1)
    else:
        print(args.root)
        ancestor = unrooted_tree.get_common_ancestor(args.root)
        unrooted_tree.set_outgroup(ancestor)
        unrooted_tree.write(format=1, outfile=args.output + "rooted_tree.nwk")
        input_tree = Tree(args.output + "rooted_tree.nwk", format = 1)
else:
    input_tree = Tree(args.tree, format = 1)

print(input_tree)

##Get list of species from input species tree & label internal nodes
tree_sp_list = []
edge = 0
sp_num = {}
sp_name_num = {}
for node in input_tree.traverse("postorder"):    
    if node.is_leaf():
        node_name = node.name
        node_num = "%d" %edge
        sp_num[node_name] = node_num
        sp_name_num[node_num] = node_name
        tree_sp_list.append(node_name)

    node.name = "%d" %edge
    edge+=1

input_tree.write(format = 1, outfile = args.output + "species_tree_label.nwk")

def tree_rates(SCO):
    '''
    Function to find fast evolving
    genes in a specified lineage
    '''
    OG_name = SCO.split('/')[-1].split('.')[0]
    SCO_tree = glob.glob(args.input + 'Gene_Trees/' + OG_name + '_tree.txt')
    SCO_tree = SCO_tree[0]
    SCO_tree_rate = check_output('phykit evolutionary_rate ' + SCO_tree, shell=True)

    SCO_rates_dist = []
    SCO_tree_rate = SCO_tree_rate.decode("utf-8")
    SCO_tree_rate = SCO_tree_rate.strip()
    SCO_rates_dist.append(SCO_tree_rate)
    
    tip_rates = check_output('phykit terminal_branch_stats -v ' + SCO_tree, shell=True)
    tip_rates = tip_rates.decode("utf-8").split('\n')
    tip_rates = tip_rates[:-1]
    tip_names = check_output('phykit tip_labels ' + SCO_tree, shell=True)
    tip_names = tip_names.decode("utf-8")
    
    sp_order = []
    for gene_name in tip_names.split('\n'):
        sp_gene_name = gene_name.split('_')[0]
        sp_order.append(sp_gene_name)
    sp_order = sp_order[:-1]
    
    species_index = sp_order.index(args.species)
    species_tip_length = tip_rates[species_index]
    fast_tip_count = 0
    for sp_rate in tip_rates:
        if sp_rate > SCO_tree_rate:
            fast_tip_count+=1
    
    if (species_tip_length > SCO_tree_rate) and (fast_tip_count <= 4):
        unix('cp ' + SCO_tree + ' ' + args.output + args.species + '_gene_stats/fast_evolving_genes/', shell=True)
            

##Get info of given species if the flag is given
def sp_stats(species):
    os.makedirs(args.output + args.species + '_gene_stats/', exist_ok=True)
    os.makedirs(args.output + args.species + '_gene_stats/fast_evolving_genes/', exist_ok=True)
    os.makedirs(args.output + args.species + '_gene_stats/multi_copy_genes/', exist_ok=True)    
    os.makedirs(args.output + args.species + '_gene_stats/species_specific_genes/', exist_ok=True)
    os.makedirs(args.output + args.species + '_gene_stats/lost_genes/', exist_ok=True)
    
    #Measure rates of species in single copy gene trees
    SCO_list = []
    for SCO in glob.glob(args.input + 'Single_Copy_Orthologue_Sequences/*fa'):
        SCO_list.append(SCO)

    if args.threads:
        Parallel(n_jobs=threads)(delayed(tree_rates)(SCO_tree) for SCO_tree in SCO_list)
    else:
        Parallel(n_jobs=1)(delayed(tree_rates)(SCO_tree) for SCO_tree in SCO_list)
        
    #Output OGs that are specific to species / duplicated just in species / or lost just in species
    taxa_count = len(sp_name_num)
    sp_multi_copy_genes = []
    sp_specific_genes = []
    sp_lost_genes = []
    for OG in glob.glob(args.input + 'Orthogroup_Sequences/*fa'):
        OG_name = OG.split('/')[-1].split('.')[0] 
        sp_list = []
        with open(OG) as f:
            for record in SeqIO.parse(f, 'fasta'):
                ID = record.id
                sp = ID.split('_')[0]##This needs to be updated to account for different header names
                sp_list.append(sp)
                
        taxa_prop = (len(set(sp_list))/taxa_count)*100
        OG_taxa_count = Counter(sp_list)
        counts = Counter(OG_taxa_count.values())
        
        if (args.species in sp_list) and (taxa_prop > 50) and (OG_taxa_count[args.species] > 1) and (len(counts) ==2) and (counts['2'] == 1):
            unix('cp ' + OG + ' ' + args.output + args.species + '_gene_stats/multi_copy_genes/', shell=True)
            sp_multi_copy_genes.append(OG_name)

        if (args.species in sp_list) and (len(set(sp_list)) == 1):
            sp_specific_genes.append(OG_name)
            unix('cp ' + OG + ' ' + args.output + args.species + '_gene_stats/species_specific_genes/', shell=True)

        if (args.species not in sp_list) and (len(set(sp_list)) == taxa_count -1):
            sp_lost_genes.append(OG_name)
            unix('cp ' + OG + ' ' + args.output + args.species + '_gene_stats/lost_genes/', shell=True)
        
    print(len(set(sp_multi_copy_genes)), 'genes in multi-copy in', args.species, 'but single in all others')
    print(len(set(sp_specific_genes)), 'genes specific to', args.species)
    print(len(set(sp_lost_genes)), 'genes lost in', args.species)
    
    
##Parse Orthgroups files
OG_origins = defaultdict(list)
single_sp_OG = 0
OG_list = []
with open(args.input + 'Orthogroups/Orthogroups.GeneCount.tsv') as f:
    species_gene_counts = defaultdict(list)
    first_line = f.readline()
    species_order = first_line.split('\t')[1:-1]
    #Check to make sure tree contains all species from orthofinder output
    if len([element for element in species_order if element not in tree_sp_list]) > 0:
        print('ERROR: Check input tree - ', list(element for element in species_order if element not in tree_sp_list), 'missing from species tree')
        sys.exit()

    sp_count = len(species_order)
    print('\nParsing output for', sp_count, 'species...\n')
    next(f)
    for line in f:
        lines = line.split('\t')
        OG = lines[0]
        OG_list.append(OG)
        sp_content = lines[1:-1]
        sp_content_counter = 0
        for gene in sp_content:
            sp_match = species_order[sp_content_counter]
            sp_content_counter+=1
            if int(gene) > 0:
                species_gene_counts[OG].append(sp_match)

    for OGs, sp_lists in species_gene_counts.items():
        sp_num_lists = []
        for sp_name in sp_lists:
            sp_node_num = sp_num[sp_name]
            sp_num_lists.append(sp_node_num)

        if len(set(sp_lists)) > 2:
            ancestor = input_tree.get_common_ancestor(sp_num_lists)
            for nodes in ancestor.traverse("preorder"):
                origin_node = nodes.name
                break
            OG_origins[origin_node].append(OGs)
        else:
            single_sp_OG+=1
            sp_node = sp_lists[0]
            OG_origins[sp_node].append(OGs)

node_counts = {}
for node, gene_list in OG_origins.items():
    try:
        sp_node_num = sp_num[node]
        gain_count = len(gene_list)
        node_counts[sp_node_num] = gain_count
    except:
        gain_count = len(gene_list)
        node_counts[node] = gain_count

##Get gene numbers per species
sp_gene_counts = {}
with open(args.input + "Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv") as f:
    head = [next(f) for x in range(2)]
    sp_list = head[0].split('\t')[1:]
    gene_nums = head[1].split('\t')[1:]
    count = 0
    for sp in sp_list:
        sp = sp.strip()
        node_num_sp = sp_num[sp]
        gene_count = gene_nums[count]
        gene_count = gene_count.strip()
        gene_count = int(gene_count)
        sp_gene_counts[node_num_sp] = gene_count
        count+=1
    
    
##Plot tree with gains
tre = toytree.tree(args.output + "species_tree_label.nwk")
for node in tre.treenode.traverse():
    node_name = node.name
    if node_name in node_counts.keys():
        values = node_counts[node_name]
        node.add_feature("genegains", values)
    else:
        node.add_feature("genegains", 0)

sizes = tre.get_node_values("genegains", True, True)
with np.errstate(divide='ignore'):
    log_sizes = np.log10(sizes)
log_10 = []
for logs in log_sizes:
    log10 = logs*10
    log_10.append(log10)


modnames = [sp_name_num[tip] for tip in tre.get_tip_labels()]
gene_counts = [sp_gene_counts[tip] for tip in tre.get_tip_labels()]
Ntips = tre.ntips
if Ntips < 30:
    canvas = toyplot.Canvas(width=800, height=750)
    ax0 = canvas.cartesian(bounds=(50, 600, 10, 700), padding=15, ymin=0, ymax=Ntips)
    ax1 = canvas.cartesian(bounds=(620, 750, 10, 700), padding=15, ymin=0, ymax=Ntips)
    tre.draw(axes=ax0, tip_labels=(modnames), node_labels=("genegains", 1, 1), node_sizes=log_10, node_colors="#99d8c9", node_style={"stroke": "black"}, tip_labels_align=True, scalebar=True, tip_labels_style={"font-size":"15px"})
    ax1.bars(np.arange(Ntips), gene_counts, along='y', color = "#bdbdbd", style={"stroke": "white", "stroke-width": 5})
    ax1.show = True
    ax1.y.show = False
    ax1.x.ticks.show = True
else:
    #canvas, axes, mark = tre.draw(width=1000, height=1200, tip_labels=(modnames), node_labels=("genegains", 1, 1), node_sizes=log_10, node_colors="#99d8c9", node_style={"stroke": "black"}, tip_labels_align=True, scalebar=True)
    canvas = toyplot.Canvas(width=1000, height=1200)
    ax0 = canvas.cartesian(bounds=(50, 800, 10, 1150), padding=15, ymin=0, ymax=Ntips)
    ax1 = canvas.cartesian(bounds=(820, 950, 10, 1150), padding=15, ymin=0, ymax=Ntips)
    tre.draw(axes=ax0, tip_labels=(modnames), node_labels=("genegains", 1, 1), node_sizes=log_10, node_colors="#99d8c9", node_style={"stroke": "black"}, tip_labels_align=True, scalebar=True, tip_labels_style={"font-size":"15px"})
    ax1.bars(np.arange(Ntips), gene_counts, along='y', color = "#bdbdbd")
    ax1.show = True
    ax1.y.show = False
    ax1.x.ticks.show = True


toyplot.pdf.render(canvas, args.output + "gene_gains_sp_tree.pdf")
toyplot.svg.render(canvas, args.output + "gene_gains_sp_tree.svg")

if args.species:
    print('\nParsing stats for', args.species)
    print('This may take some time...\n')
    sp_stats(args.species)
                
print('\nSuccess! Output found in', args.output)
