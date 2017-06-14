#!/usr/bin/env python

import argparse
#from Bio import Phylo

def find_all(string, pattern):
    return [i for i in range(len(string)) if string.startswith(pattern, i)]

def filter_num(ref_file, threshold):
    with open(ref_file) as f:
        f.readline()
        for line in f:
            line = line.strip()
            oc = find_all(line, 'BR')
            if len(oc) >= threshold:
                print line
 
def get_children(tree):
    parents = {}
    for clade in tree.find_clades(order='level'):
        parents[clade.name] = map(lambda x: x.name, clade.get_terminals())
    return parents

# def search_tree(ref_file, tree, ancs):
#     parents = get_children(tree)
#     with open(ref_file) as f:
#         species = f.readline().strip().split()[1:]
#         #print species
#         for line in f:
#             line = line.strip()
#             data = line.split()
#             brs = [i for i,val in enumerate(data) if val=='BR']
#             #print 'brs:', brs
#             oc = map(lambda x:species[x],brs)
#             #print 'oc:', oc
#             if oc:
#                 com_anc = tree.common_ancestor(oc)
#                 if com_anc.name in ancs and len(parents[com_anc.name]) == len(oc):
#                     print com_anc.name+':', line
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='breakpoints table')
    parser.add_argument('--filter_num', help='filter rows containg number of breakpoints')
    # parser.add_argument('--search_tree', help='filter newick tree')
    # parser.add_argument('--key_ancs', nargs='+', help='list of ancestral genomes')
    args = parser.parse_args()

    if args.filter_num:
        filter_num(args.file, int(args.filter_num))

    # if args.search_tree:
    #     tree = Phylo.read(args.search_tree, 'newick')
    #     search_tree(args.file, tree, args.key_ancs)

