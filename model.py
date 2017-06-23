from collections import Counter


class Chromosome:
    def __init__(self, seq_id, size, description):
        self.seq_id = seq_id
        self.size = size
        self.description = description

    def get_specie(self):
        return self.description.split('.')[0]

    def print_out(self):
        print self.seq_id, self.size, self.description


class Entry:
    def __init__(self, seq_id, strand, start, end, length):
        self.seq_id = seq_id
        self.strand = strand
        self.block_id = -1
        if start > end:
            self.start = end
            self.end = start
        else:
            self.start = start
            self.end = end
        self.length = length

    def get_seq_id(self):
        return self.seq_id

    def get_start(self):
        return self.start

    def get_end(self):
        return self.end

    def set_block_id(self, block_id):
        self.block_id = block_id

    def get_block_id(self):
        return self.block_id

    def get_specie(self):
        return self.seq_id.split('.')[0]

    def get_chrom(self):
        return self.seq_id.split('.')[1]

    def equals(self, e):
        return self.seq_id == e.seq_id and self.strand == e.strand\
                and self.start == e.start and self.end == e.end\
                and self.length == e.length

    def equals_to_list(self, e_list):
        for x in e_list:
            if not self.equals(x):
                return False
        return True

    def print_out(self):
        print 'seq_id:', self.seq_id, 'block_id:', self.block_id,
        'strand:', self.strand, 'start:', self.start, 'end:', self.end,
        'length:', self.length

    def print_out_short(self):
        print ' '.join(map(str, [self.get_seq_id(), self.start, self.end]))


class Block:
    def __init__(self, id, entries):
        self.id = id
        self.entries = entries
        for e in self.entries:
            e.set_block_id(id)

    def print_out(self):
        print self.id
        for e in self.entries:
            e.print_out()

    def get_species(self):
        return set(map(lambda x: x.get_specie(), self.entries))


class BED_Entry:
    def __init__(self, genome, chrom, start, end):
        self.genome = genome
        self.chrom = chrom
        self.start = start
        self.end = end

    def to_string(self):
        return '_'.join(map(str, [self.genome + '.' + self.chrom,
                                  self.start, self.end]))

    def print_out(self):
        print ' '.join(map(str, [self.genome + '.' + self.chrom,
                                 self.start, self.end]))


def parse_chromosomes(f):
    chroms = {}
    with open(f) as blocks_file:
        blocks_file.readline()
        for line in blocks_file:
            line = line.strip().split()
            if len(line) != 3:
                return chroms
            chroms[int(line[0])] = Chromosome(int(line[0]), int(line[1]),
                                              line[2])
    return chroms


SPLITTER = '-----------------------------------'


def parse_blocks(f, count_c=False, skip_dups=False):
    count_chrs = {}
    max_block_id = 0
    with open(f) as blocks_file:
        blocks_section = False
        entries = []
        blocks = []
        id = ''
        for line in blocks_file:
            line = line.strip()
            if 'Block #' in line:
                blocks_section = True
            if blocks_section:
                if 'Block #' in line:
                    if id:
                        c = Counter(map(lambda x: x.get_specie(), entries))
                        if not skip_dups or not filter(lambda x: x > 1,
                                                       c.values()):
                            blocks.append(Block(id, entries))
                    id = int(line[7:])
                    if id > max_block_id:
                        max_block_id = id
                    entries = []
                    continue
                if 'Seq_id' in line or SPLITTER in line:
                    continue
                line = line.split()
                seq_id = line[0]
                # if seq_id.find('random') > -1 or seq_id.find('hap') > -1:
                #    continue
                entries.append(Entry(seq_id, line[1], int(line[2]),
                                     int(line[3]), int(line[4])))
                if count_c:
                    if seq_id in count_chrs.keys():
                        count_chrs[seq_id] += 1
                    else:
                        count_chrs[seq_id] = 1
        c = Counter(map(lambda x: x.get_specie(), entries))
        if not skip_dups or not filter(lambda x: x > 1, c.values()):
            blocks.append(Block(id, entries))
    if count_c:
        return blocks, count_chrs, max_block_id
    else:
        return blocks
