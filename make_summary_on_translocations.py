#!/usr/bin/env python

import sys
import os
import itertools
import argparse

class TranslocationEntry :
    related_chr1 = None
    related_chr2 = None
    specie = None
    blocks = []

    def __init__ (self, related_chr1, related_chr2, specie, blocks):
        self.related_chr1 = related_chr1
        self.related_chr2 = related_chr2
        self.specie = specie
        self.blocks = blocks

    def add_block(self, x):
        self.blocks.append(x)


def overlap_common_entries(entries, type) :
    all_blocks = map(lambda x: x.blocks, entries)
    all_blocks = set(itertools.chain(*all_blocks))
    entries_by_key_all = []
    merged = []
    #make list of lists of entries that share one common key
    for key in all_blocks:
        entries_by_key = filter(lambda x: key in x.blocks, entries)
        entries_by_key_all.append(entries_by_key)
    #merge list that are contained by larger ones
    while (entries_by_key_all) :
        key_list = entries_by_key_all[-1]
        related = filter(lambda x: set(x) <= set(key_list) or set(x) >= set(key_list), entries_by_key_all)
        for x in related:
             entries_by_key_all.remove(x)
        related = set(itertools.chain(*related))
        merged.append(related)

    i = 0
    for x in merged:
        if len(x) < 2:
            continue
        i += 1
        print type + ' no.', i
        for y in x:
            print y.specie, y.related_chr1, y.related_chr2, y.blocks
        print


def parse_translocation(dir_name, type):
    summary = []
    for name in os.listdir(dir_name):
        file_name = os.path.join(dir_name, name)
        with open(file_name) as f:
            specie = f.readline().strip()
            chr_name1 = None
            chr_name2 = None
            related_blocks = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if type in line:
                    if chr_name1 and chr_name2:
                        if not chr_name1:
                            chr_name1 = chr_name2
                        else:
                            chr_name1, chr_name2 = sorted([chr_name1, chr_name2])
                        summary.append(TranslocationEntry(chr_name1, chr_name2, specie, related_blocks))
                        related_blocks = []
                    chr_name1 = line.split("from chromosome ")
                    if len(chr_name1) > 1:
                        chr_name1 = chr_name1[1]
                    else:
                        chr_name1 = None
                elif 'seq_id: ' in line:
                    data = line.split('seq_id: ')[1]
                    data = data.split(' block_id: ')
                    chr_name2 = data[0]
                    b = data[1]
                    b = b.split(' strand: ')[0]
                    related_blocks.append(int(b))
                    chr_name2 = chr_name2.split('.')[1]
    return summary

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', required=True, help='file with rearrangements')
    parser.add_argument('--type', required=True, help='translocation/transposition')
    args = parser.parse_args()
    summary = parse_translocation(args.dir, args.type)
    overlap_common_entries(summary, args.type)
