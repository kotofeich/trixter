import model
import itertools


def intersect(entry, bed_entries):
    intersected_bed_entries = []
    for b in bed_entries:
        if (entry.get_specie() == b.genome or b.genome == '') and entry.get_chrom() == b.chrom:
            if not (b.end < entry.start or entry.end < b.start):
                intersected_bed_entries.append(b)
    return intersected_bed_entries


'''
traverses blocks and collects all the entries
related to the specie
'''
def get_specie_entries(blocks, specie):
    specie_entries = []
    for b in blocks:
        # i = 0
        for e in b.entries:
            if e.get_specie() == specie:
                specie_entries.append(e)
    return specie_entries

def reorder_specie(specie1, entries2, name1, name2):
    specie2_grouped = []
    # group entries in specie2 according to the order of blocks on chromosomes in specie1
    visited_blocks = set()
    order = []
    for sp in specie1:
        specie2_grouped.append([])
        # unpaired entries are erased
        unpaired_entries = set([])
        for y in sp:
            if not order:
                order.append(y.seq_id)
            if order and order[-1] != y.seq_id:
                order.append(y.seq_id)
            if y.block_id in visited_blocks:
                continue
            if len(filter(lambda a: a.block_id == y.block_id, sp)) > 1:
                raise Exception('duplicated block', y.block_id, 'in', name1)
            c = filter(lambda x: x.block_id == y.block_id, entries2)
            visited_blocks.add(y.block_id)
            if len(c) == 1:
                specie2_grouped[-1].append(c)
            elif len(c) > 1:
                raise Exception('duplicated block', y.block_id, 'in', name1)
            elif not c:
                unpaired_entries.add(y)
                print 'no such blocks ', y.block_id, 'in specie', name2
                # now let's order the blocks that are duplicated on the same chromosome
        for y in unpaired_entries:
            sp.remove(y)
    specie2_rear = []
    for e in specie2_grouped:
        specie2_rear.append(list(itertools.chain(*e)))
    return specie2_rear

def thread_specie_genome(specie_entries):
    genome = []
    chromosomes_names = set(map(lambda x: x.get_seq_id(), specie_entries))
    chromosomes = []
    # group entries by chromosomes
    for c in chromosomes_names:
        c_entries = filter(lambda x: x.get_seq_id() == c, specie_entries)
        chromosomes.append(c_entries)
    # sort entries in chromosomes by position
    sorted_chromosomes = []
    for c in chromosomes:
        sorted_c = sorted(c, key=lambda x: x.get_start())
        sorted_chromosomes.append(sorted_c)
    return sorted_chromosomes


def print_out_genome_thread(entries, path=None):
    i = 0
    if path:
        f = open(path, 'w')
    for c in entries:
        i += 1
        if path:
            f.write(str(i) + '\n')
        else:
            print i
        for e in c:
            s = 'seq_id: ' + str(e.get_seq_id()) + ' block_id: ' + str(e.get_block_id()) + ' strand: ' \
                + str(e.strand) + ' start: ' + str(e.get_start()) + ' end: ' + str(e.get_end())
            if path:
                f.write(s + '\n')
            else:
                print s
    if path:
        f.close()


def group_by_ref(ref_genome, target_genome):
    target_genome_grouped = []
    visited_blocks = set()
    for sp in ref_genome:
        target_genome_grouped.append([])
        # unpaired entries are erased
        unpaired_entries = set([])
        for y in sp:
            if y.get_block_id() in visited_blocks:
                continue
            c = filter(lambda x: x.get_block_id() == y.get_block_id(), list(itertools.chain(*target_genome)))
            visited_blocks.add(y.get_block_id())
            if len(c) == 1:
                target_genome_grouped[-1].append(c[0])
            elif not c:
                unpaired_entries.add(y)
            else:
                raise Exception('duplicated block!')
            for y in unpaired_entries:
                sp.remove(y)
    return target_genome_grouped


def get_neighbors(c, e):
    ind = c.index(e)
    if ind:
        this_prev_blocks_id = c[ind - 1].block_id
    else:
        this_prev_blocks_id = None
    if ind < len(c) - 1:
        this_next_blocks_id = c[ind + 1].block_id
    else:
        this_next_blocks_id = None
    return (this_prev_blocks_id, this_next_blocks_id)


# normalization means we revert all the negative-strand blocks of the chromosome in specie1
# and change the strand of the corresponding block in specie2
# this is needed in order to search for reversals only in specie2 related to specie1
def normalize(specie1, specie2):
    for i in range(len(specie1)):
        c1 = specie1[i]
        c2 = specie2[i]
        if not c1 or not c2:
            print 'skipping empty chromosome'
            continue
        for j in range(len(c1)):
            if c1[j].strand == '-':
                c1[j].strand = '+'
                if c2[j].strand == '-':
                    c2[j].strand = '+'
                elif c2[j].strand == '+':
                    c2[j].strand = '-'
        specie1[i] = c1
        specie2[i] = c2
    return specie1, specie2


def filter_unsplitted_chromosomes(blocks, count_chrs, sps):
    upd_blocks = []
    for b in blocks:
        entries = b.entries
        upd_entries = []
        upd_species = set()
        for e in entries:
            specie = e.get_specie()
            if specie in sps:
                if count_chrs[e.seq_id] > 1:
                    upd_entries.append(e)
                    upd_species.add(specie)
        # also count duplications?
        # if so then only blocks when both chromosomes are split counted
        if len(upd_entries) >= len(sps) and len(upd_species) == len(sps):
            # if so then counted also those blocks that partly split but in some
            # species it can be the whole scaffold
            # if upd_entries:
            upd_blocks.append(model.Block(b.id, upd_entries))
    return upd_blocks


def filter_absent_species(blocks, sps):
    upd_blocks = []
    for b in blocks:
        entries = b.entries
        upd_species = []
        for e in entries:
            specie = e.get_specie()
            if specie in sps:
                upd_species.append(specie)
        if len(upd_species) == len(sps):
            upd_blocks.append(b)
    return upd_blocks


# must be fixed: in case of duplications in a genome
# there can be ambiguities in prev entries
def find_prev_block_in_specie(entry, specie):
    for c in specie:
        find = filter(lambda x: x.block_id == entry.block_id and x.start == entry.start, c)
        if len(find) > 1:
            raise Exception('duplicated entry!')
        if find:
            l = c.index(find[0])
            if l == 0:
                return None
            return c[l - 1]
    raise Exception('No such block! ', entry.block_id)


# must be fixed: in case of duplications in a genome
# there can be ambiguities in next entries
def find_next_block_in_specie(entry, specie):
    for c in specie:
        find = filter(lambda x: x.block_id == entry.block_id and x.start == entry.start, c)
        if len(find) > 1:
            raise Exception('duplicated entry!')
        if find:
            l = c.index(find[0])
            if l == len(c) - 1:
                return None
            return c[l + 1]
    raise Exception('No such block! ', entry.block_id)
