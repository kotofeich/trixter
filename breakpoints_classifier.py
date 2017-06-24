#!/usr/bin/env python
# -*- coding: utf-8 -*-

import utils
import itertools
from collections import Counter, defaultdict


class Decision:
    true = True
    false = False
    noinfo = None


def get_set_entries(blocks):
    list_of_list_entries = map(lambda x: x.entries, blocks)
    set_of_entries = set(list(itertools.chain(*list_of_list_entries)))
    set_of_species = set(map(lambda x: x.get_specie(), set_of_entries))
    return set_of_species


def get_blocks_ids(blocks):
    return set(map(lambda x: x.id, blocks))


def create_indices(species, threaded_genomes):
    ind = defaultdict(list)
    for g in threaded_genomes.keys():
        chrs = threaded_genomes[g]
        for k in range(len(chrs)):
            for i in range(len(chrs[k])):
                prev_id = -1 if (i == 0) else chrs[k][i-1].block_id
                next_id = -2 if (i == len(chrs[k]) - 1) else \
                    chrs[k][i+1].block_id
                ind[(chrs[k][i].block_id, g)].append((prev_id, next_id))
    return ind


# doesnt work for duplications (will rewrite the previous entry)
def build_features_index(blocks):
    features = {}
    for b in blocks:
        if len(set(map(lambda x: x.get_specie(), b.entries))) > \
           len(map(lambda x: x.get_specie(), b.entries)):
            print
            raise Exception('cant handle duplications in ' +
                            build_features_index.__name__)
        for e in b.entries:
            features[(e.get_specie(), e.block_id)] = [e.get_chrom(), e.start,
                                                      e.end]
    return features


def get_features(specie, block_id, features):
    if block_id == '-3':
        return ' '.join([specie, 'duplication'])
    return ' '.join(map(str, features[(specie, block_id)]))


def build_neighbours(block_inds, index):
    # just linearize two-dimensional data
    # [(prev1, next1), (prev2, next2)] -> [prev1, next1, prev2, next2]
    entries_num = 0
    dupls_num = 0
    neighbours = []
    species_status = {}
    for ind in block_inds:
        entries_num += len(index[ind])
        if len(index[ind]) > 1:
            dupls_num += len(index[ind])
            index[ind] = [(-3, -3)]
            species_status[ind[1]] = 'DUP'
        # beware of dupl!
        neighb = index[ind][0]
        neighbours.append(neighb[0])
        neighbours.append(neighb[1])
    return dupls_num, entries_num, neighbours, species_status


def process_at_index(index, ind, allowable, features, species_status, nodef,
                     print_table):
    # beware of dupl!
    prev, next = index[ind][0]
    if prev in nodef and next in nodef:
        # the whole block is a full scaffold in the specie
        species_status[ind[1]] = 'END'
        return
    if prev in allowable and next in allowable:
        species_status[ind[1]] = '-'
        return
    if prev not in allowable:
        species_status[ind[1]] = 'BR'
        if not print_table:
            print 'breakpoint', ind[1], prev,
            get_features(ind[1], prev, features), '-', ind[0],
            get_features(ind[1], ind[0], features)
    if next not in allowable:
        species_status[ind[1]] = 'BR'
        if not print_table:
            print ind[1]
            print ind[0]
            print next
            print 'breakpoint', ind[1], ind[0],
            get_features(ind[1], ind[0], features), '-', next,
            get_features(ind[1], next, features)


def process_block_neighborhood(neighbours, block_inds, species_status,
                               print_table, features, index):
    # sort by popularity in descending order
    # and leave only non-ending
    c = Counter(neighbours).most_common()
    c = filter(lambda x: x[0] != -2 and x[0] != -1 and x[0] != -3, c)
    # if len is less or equal than two (most popular from left and from right),
    # then breakpoint is likely to be caused by assembly incompleteness
    if len(c) > 2:
        if c[2][1] == c[1][1]:
            if not print_table:
                print 'cant distinguish two most common!'
                print
            for ind in block_inds:
                if not ind[1] in species_status.keys():
                    # could not resolve breakpoint
                    species_status[ind[1]] = 'NA'
            return
        first_common = c[0][0]
        second_common = c[1][0]
        nodef = set([-1, -2])
        allowable = set([-1, -2, first_common, second_common])
        for ind in block_inds:
            process_at_index(index, ind, allowable, features, species_status,
                             nodef, print_table)


def run(blocks, print_table=False):
    species = sorted(list(get_set_entries(blocks)))
    threaded_genomes = {}
    for sp in species:
            entries = utils.get_specie_entries(blocks, sp)
            threaded_genomes[sp] = utils.thread_specie_genome(entries)
    index = create_indices(species, threaded_genomes)
    blocks_ids = get_blocks_ids(blocks)
    blocks_num = 0

    features = {}
    if not print_table:
        features = build_features_index(blocks)
    if print_table:
        header = '\t'.join(['breakpoint block']+species)
        print header
    for b in blocks_ids:
        blocks_num += 1
        block_inds = filter(lambda x: x[0] == b, index)
        if not print_table:
            print 'block_id:', b
            for ind in block_inds:
                print index[ind]
        dupls_num, entries_num, neighbours, species_status = build_neighbours(
            block_inds, index)
        process_block_neighborhood(neighbours, block_inds, species_status,
                                   print_table, features, index)
        if not print_table:
            print
        if print_table:
            l = str(b)
            if not species_status.keys():
                # in case breakpoint is caused by assembly incompleteness
                continue
            for e in species:
                if e not in species_status.keys():
                    l += '\t' + 'not-in-block'
                else:
                    l += '\t' + species_status[e]
            print l
    if not print_table:
        print 'STAT number of blocks:', blocks_num
        print 'STAT number of entries:', entries_num
        print 'STAT number of dupls (among entries):', dupls_num
        print 'STAT rate of duplications:', float(dupls_num)/entries_num
