from collections import Counter, defaultdict
import itertools
import math
import utils

def check_order(list_a, sorted_list):
    #join neighbours
    neighbours = defaultdict(set)
    sorted_neighbours = defaultdict(set)
    for i in range(len(list_a)):
        if i > 0:
            neighbours[list_a[i]].add(list_a[i-1])
            sorted_neighbours[sorted_list[i]].add(sorted_list[i-1])
        else:
            neighbours[list_a[i]].add(None)
            sorted_neighbours[sorted_list[i]].add(None)
        if i < len(list_a) - 1:
            neighbours[list_a[i]].add(list_a[i+1])
            sorted_neighbours[sorted_list[i]].add(sorted_list[i+1])
        else:
            neighbours[list_a[i]].add(None)
            sorted_neighbours[sorted_list[i]].add(None)

    '''
    print 'n'
    for e in neighbours:
        e.print_out()
        print '-'
        for x in neighbours[e]:
            x.print_out()
        print
    print 'sorted n'
    for e in sorted_neighbours:
        e.print_out()
        print '-'
        for x in sorted_neighbours[e]:
            x.print_out()
        print
    '''
    transpositions = []
    reported = set()
    #check if neighbours are different
    for e in neighbours:
        diff = sorted_neighbours[e] - neighbours[e]
        if diff:
            for x in diff:
                if not (e,x) in reported and not (x,e) in reported:
                    transpositions.append((e,x))
                    reported.add((e,x))
    return transpositions

'''
we need it for two reasons:
1. in order to count rearrangements we look at the previous entry and
check if it's rearranged. sometimes we can find multiple possible previous entries.
in this case we need to know which one is previous in this case
2. in order to be able to report the breakpoints
'''
def get_previous_entries(list_entries, c):
    rearrangement_ids = []
    for e in list_entries:
        rearrangement_ids.append(c.index(e))
    rearrangement_ids = map(lambda x: x-1, rearrangement_ids)
    rearrangement_prev = []
    for i in rearrangement_ids:
        if i == -1:
            rearrangement_prev.append(None)
        else:
            rearrangement_prev.append(c[i])
    return rearrangement_prev

def get_next_entries(list_entries, c):
    rearrangement_ids = []
    for e in list_entries:
        rearrangement_ids.append(c.index(e))
    rearrangement_ids = map(lambda x: x+1, rearrangement_ids)
    rearrangement_prev = []
    for i in rearrangement_ids:
        if i == len(c):
            rearrangement_prev.append(None)
        else:
            rearrangement_prev.append(c[i])
    return rearrangement_prev

'''
transposition is a change of the genomic location on the same chromosome
idea: we report a couple of entries (e1,e2) that are not neighbours in sorting
the entries according to the order of the other genome
'''
def check_transpositions(c):
    if not c:
        return ([],[],[])
    transpositions = []
    seq_ids = map(lambda x: x.get_seq_id(), c)
    seq_ids = set(seq_ids)
    #consider blocks separately for each chromosome
    for s in seq_ids:
        c_seq_id = filter(lambda x: x.get_seq_id() == s, c)
        c_seq_id_sorted = list(c_seq_id)
        c_seq_id_sorted.sort(key=lambda x: x.get_start())
        if (c_seq_id == c_seq_id_sorted) or (c_seq_id == c_seq_id_sorted[::-1]):
            continue
        '''
        print 'aj according fc:'
        for x in c_seq_id:
            x.print_out() 
        print 'aj sorted regularly:'
        for x in c_seq_id_sorted:
            x.print_out()
        '''
        transpositions += check_order(c_seq_id, c_seq_id_sorted)
    return transpositions

'''
translocation is a change of the genomic position to another chromosome
the main chromosome is the one that corresponds to less translocations (by length)
generally translocations can be counted without information about previous entries
because they are already grouped into lists of unbroken segments of blocks.
however we report translocations in the same manner as reversals and transpositions
for uniformity and for being able to report the breakpoints if needed
'''
def check_translocations(c):
    if not c:
        return [], []
    c_seq_ids =[]
    c_sorted = sorted(c, key=lambda x:x.get_seq_id())
    for e in itertools.groupby(c_sorted, lambda x: x.get_seq_id()):
        c_seq_ids.append(list(e[1]))
    lengths = map(lambda y: sum(map(lambda x: math.fabs(int(x.get_end())-int(x.get_start())), y)), c_seq_ids)
    ls = list(zip(lengths, c_seq_ids))
    ls_sorted = sorted(ls, key=lambda x: x[0])
    translocations = map(lambda x: x[1], ls_sorted[:-1])
    #main_chrom = (ls_sorted[-1][1][0].get_chrom(), ls_sorted[-1][1][0].block_id)
    main_chrom = ls_sorted[-1][1][0].get_chrom()
    '''
    print 'start'
    for e in ls_sorted:
        print '----'
        for x in e[1]:
            x.print_out()
        print '----'
        '''
    return main_chrom, translocations

'''
reversal is the change of strand
here we use the normalization of two genomes. i.e. we brought the genome of specie1
to the form where all the entries have '+' strand, changing also the strand of
the corresponding entries in specie2
every '-' is called reversal
but if the whole 'chromosome' is '-' than nothing is reversed
'''
def check_reversals(c):
    if not c:
        return []
    c_rev = filter(lambda x: x.strand == '-', c)
    if len(c_rev) == len(c) :
        return []
    if c_rev:
        rev_prev = get_previous_entries(c_rev,c)
        return zip(rev_prev,c_rev)
    return []
'''
If block appears in genome several times then it's a duplication
returns species entries that belong to duplicated blocks
'''
def check_duplications(c, blocks, specie):
    l = utils.get_specie_entries(blocks,specie)
    cnt = Counter(map(lambda x: x.block_id, l))
    dups = filter(lambda x: cnt[x.block_id] > 1, c)
    if dups:
        dup_prev = get_previous_entries(dups,c)
        return zip(dup_prev,dups)
    else:
        return []
