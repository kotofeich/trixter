#!/usr/bin/env python

import os
import argparse
import itertools
from collections import defaultdict

import model
import utils
from blocks_to_paths_processor import BlocksToPathsProcessor
import rearrangements_type
import breakpoints_classifier

def print_out_genome_thread(entries):
    #with open(file_name,'w') as f:
        i = 0
        for c in entries:
            i += 1
            #f.write(str(i)+'\n')
            print i
            for e in c:
                 #f.write('seq_id: ' + str(e.seq_id) + ' block_id: ' + str(e.block_id) + ' strand: '\
                 #+ str(e.strand) + ' start: ' + str(e.start) + ' end: ' + str(e.end) + '\n')
                 print 'seq_id: ' + str(e.seq_id) + ' block_id: ' + str(e.block_id) + ' strand: '\
                 + str(e.strand) + ' start: ' + str(e.start) + ' end: ' + str(e.end)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='blocks_coords.txt')
    parser.add_argument('--report_transpositions', action='store_true', help='report transpositions in specie2 related to specie1')
    parser.add_argument('--report_translocations', action='store_true', help='report translocations in specie2 related to specie1')
    parser.add_argument('--report_reversals', action='store_true', help='report reversals in specie2 related to specie1')
    parser.add_argument('--species', nargs='+', help='species to check')
    parser.add_argument('--rename_duplications', action='store_true')
    
    parser.add_argument('--classify_breakpoints', action='store_true', help='find out which species contain breakpoint')
    parser.add_argument('--print_table', action='store_true', help='not reporting themself but the list of species that contain it')

    parser.add_argument('--print_out_genomes', action='store_true', help='prints out genomes of --species in terms of blocks')
    parser.add_argument('--filter', help='filter blocks for regions mentioned in bed file')

    args = parser.parse_args()
    chroms = model.parse_chromosomes(args.file)
    blocks, count_chrs, max_block_id = model.parse_blocks(args.file, True)

    if args.filter:
        f_blocks = utils.filter_bed(blocks, args.filter)
    elif args.classify_breakpoints:
        if args.print_table:
            breakpoints_classifier.run(blocks, True)
        else:
            breakpoints_classifier.run(blocks, False)

    elif args.report_translocations or args.report_transpositions or args.report_reversals:
            if len(args.species) != 2:
                raise Exception("Can evaluate rearrangements only between two species")
            blocks = utils.filter_unsplitted_chromosomes(blocks, count_chrs, args.species)
            #get all the entries from specie1
            entries = utils.get_specie_entries(blocks, args.species[0])
            #sort entries by chromosomes for specie1
            specie1 = utils.thread_specie_genome(entries)
            #Tfor testing purposes
            #Ttest_path = '/Users/admin/projects/felidae_comp/analysis/synteny/solenodon/1000/'
            #Tprint_out_genome_thread(args.species[0],specie1,os.path.join(test_path,'tmp1'))
            entries2 = utils.get_specie_entries(blocks, args.species[1])
            #Tprint_out_genome_thread(args.species[1],utils.thread_specie_genome(entries),os.path.join(test_path,'tmp2'))
            #Texit()
            specie2 = utils.thread_specie_genome(entries2)
            if args.rename_duplications:
                specie1, renamed_prev_species, min_next_block_id = utils.rename_duplications(specie1, [], max_block_id+1)
                specie2, renamed_prev_species, min_next_block_id = utils.rename_duplications(specie2, renamed_prev_species, min_next_block_id)
            specie2_grouped = []
            #group entries in specie2 according to the order of blocks on chromosomes in specie1
            visited_blocks = set()
            for sp in specie1:
                specie2_grouped.append([])
                for y in sp:
                    if y.block_id in visited_blocks:
                        continue
                    if len(filter(lambda a:a.block_id==y.block_id, sp)) > 1:
                        raise Exception('duplicated block in', args.species[0])
                    c = filter(lambda x: x.block_id == y.block_id, entries2)
                    visited_blocks.add(y.block_id)
                    if len(c) == 1:
                        specie2_grouped[-1].append(c)
                    elif len(c) > 1:
                        raise Exception('duplicated block in', args.species[1])
                    elif not c:
                        print 'no such blocks ', y.block_id, 'in specie', args.species[1]
                    '''
                    #now let's order the blocks that are duplicated on the same chromosome
                        c_grouped_same_chrom = [list(v) for k,v in itertools.groupby(c,key=lambda x:x.seq_id)]
                        c_grouped_same_chrom = map(lambda l: sorted(l, key=lambda y: y.start), c_grouped_same_chrom)
                        c_grouped_same_chrom = list(itertools.product(*c_grouped_same_chrom))
                        specie2_grouped[-1] += c_grouped_same_chrom
                    '''
            specie2_rear = []
            cnt_empty = 0
            for e in specie2_grouped:
                p = BlocksToPathsProcessor.search_paths(e)
                if not p:
                    cnt_empty += 1
                specie2_rear.append(p)
            print 'unresolved paths (chromosomes):', cnt_empty
            specie1,specie2_rear = utils.normalize(specie1, specie2_rear)
            for c in specie2_rear:
                if args.report_transpositions:
                    trp = rearrangements_type.check_transpositions(c)
                    this_trp = []
                    to_start = -1
                    to_end = -2
                    for e in trp:
                        this_prev = e[0]
                        this_trp.append(e[1])
                        this_next = e[2]
                        #start is the block before the trasposition
                        if not this_prev in map(lambda x: x[1], trp):
                            to_start = this_prev
                        if not this_next in map(lambda x: x[1], trp):
                            to_end = this_next
                        if to_start != -1 and to_end != -2:
                            from_start = utils.find_prev_block_in_specie(this_trp[0],specie2)
                            from_end = utils.find_next_block_in_specie(this_trp[-1], specie2)
                            l = str(from_start.block_id)+'\t' if from_start else 'None\t'
                            l += str(from_end.block_id)+'\t' if from_end else 'None\t'
                            l += str(to_start.block_id)+'\t' if to_start else 'None\t'
                            l += str(to_end.block_id) if to_end else 'None'
                            for t in this_trp:
                                l += '\t'+str(t.block_id)
                            l += '\ttransposition'
                            print l
                            to_start = -1
                            to_end = -2
                            this_trp = []
                    #print 'transposition:'
                    #for t in trp:
                    #    t[1].print_out()
                if args.report_translocations:
                    main_chrom, trl = rearrangements_type.check_translocations(c)
                    for e in trl:
                        #if not this_prev in all_translocated_entries:
                        #    count_trl += 1
                        #print 'whole chromosome:'
                        #for x in c:
                        #    x.print_out()
                        print 'translocation: from chromosome', main_chrom
                        for x in e:
                            x.print_out()
                    if trl:
                        print 'overall translocations:', len(trl)
                if args.report_reversals:
                    count_rev = 0
                    rev = rearrangements_type.check_reversals(c)
                    #print 'whole chromosome'
                    #for x in c:
                    #    x.print_out()
                    for e in rev:
                        this_prev = e[0]
                        this_rev = e[1]
                        #count reversal only once if
                        #it occured in neighbouring blocks
                        if not this_prev in map(lambda x: x[1], rev):
                            count_rev += 1
                        print 'reversal:',
                        this_rev.print_out()
                    if count_rev:
                        print 'overall reversals', count_rev
    elif args.print_out_genomes :
        #blocks = utils.filter_unsplitted_chromosomes(blocks, count_chrs, args.species)
        for sp in args.species:
            entries = utils.get_specie_entries(blocks, sp)
            specie_genome = utils.thread_specie_genome(entries)
            print_out_genome_thread(specie_genome)

