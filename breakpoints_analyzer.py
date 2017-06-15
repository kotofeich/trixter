#!/usr/bin/env python

import argparse

import model
import utils
import rearrangements_type
import breakpoints_classifier


def process_transpositions(c):
    trp = rearrangements_type.check_transpositions(c)
    this_trp = []
    to_start = -1
    to_end = -2
    for e in trp:
        if not e:
            continue
        this_prev = e[0]
        this_trp.append(e[1])
        this_next = e[2]
        # start is the block before the trasposition
        if not this_prev in map(lambda x: x[1], trp):
            to_start = this_prev
        if not this_next in map(lambda x: x[1], trp):
            to_end = this_next
        if to_start != -1 and to_end != -2:
            print 'transposition'
            for t in this_trp:
                t.print_out()
            to_start = -1
            to_end = -2
            this_trp = []


def process_translocations(c):
    main_chrom, trl = rearrangements_type.check_translocations(c)
    for e in trl:
        print 'translocation: from chromosome', main_chrom
        for x in e:
            x.print_out()
    if trl:
        print 'overall translocations:', len(trl)


def process_reversals(c):
    count_rev = 0
    rev = rearrangements_type.check_reversals(c)
    for e in rev:
        if not e:
            continue
        this_prev = e[0]
        this_rev = e[1]
        # count reversal only once if
        # it occured in neighbouring blocks
        if not this_prev in map(lambda x: x[1], rev):
            count_rev += 1
        print 'reversal:',
        this_rev.print_out()
    if count_rev:
        print 'overall reversals', count_rev


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--file', required=True, help='blocks_coords.txt')
    rearrangements_group = parser.add_argument_group()
    type_group = rearrangements_group.add_mutually_exclusive_group()
    type_group.add_argument('--report_transpositions', action='store_true', help='report transpositions in specie2 related to specie1')
    type_group.add_argument('--report_translocations', action='store_true', help='report translocations in specie2 related to specie1')
    type_group.add_argument('--report_reversals', action='store_true', help='report reversals in specie2 related to specie1')
    type_group.add_argument('--report_reorganized_genome', action='store_true', help='print out genome of specie2 reordered according to the specie1')
    type_group.add_argument('--report_breakpoints', action='store_true', help='output list of breakpoints in the specie1 genome related to specie2')
    rearrangements_group.add_argument('--species', nargs=2, help='two species to evaluate rearrangements')

    breakpoints_group = parser.add_argument_group()
    breakpoints_group.add_argument('--classify_breakpoints', action='store_true', help='find out which species contain breakpoint')
    breakpoints_group.add_argument('--print_table', action='store_true', help='not reporting themself but the list of species that contain it')


    io_group = parser.add_argument_group()
    io_group.add_argument('--print_genomes', nargs='+', help='prints out specified genomes')

    args = parser.parse_args()
    chroms = model.parse_chromosomes(args.file)
    blocks, count_chrs, max_block_id = model.parse_blocks(args.file, count_c=True, skip_dups=True)

    if args.classify_breakpoints:
        if args.print_table:
            breakpoints_classifier.run(blocks, True)
        else:
            breakpoints_classifier.run(blocks, False)

    elif args.report_translocations or args.report_transpositions or \
            args.report_reversals or args.report_reorganized_genome or\
                args.report_breakpoints:
        if not args.species:
            print 'Choose species to find rearrangements --species'
            parser.print_help()
            exit()
        blocks = utils.filter_unsplitted_chromosomes(blocks, count_chrs, args.species)
        #get all the entries from specie1
        entries = utils.get_specie_entries(blocks, args.species[0])
        #sort entries by chromosomes for specie1
        specie1 = utils.thread_specie_genome(entries)
        entries2 = utils.get_specie_entries(blocks, args.species[1])
        specie2_rear = utils.reorder_specie(specie1, entries2, args.species[0], args.species[1])
        specie1,specie2_rear = utils.normalize(specie1, specie2_rear)
        if args.report_reorganized_genome:
            utils.print_out_genome_thread(specie2_rear, blocks, print_ref_id=True)
        elif args.report_breakpoints:
            utils.report_breakpoints(specie1, specie2_rear)
        else:
            for c in specie2_rear:
                if args.report_transpositions:
                    process_transpositions(c)
                if args.report_translocations:
                    process_translocations(c)
                if args.report_reversals:
                    process_reversals(c)
    elif args.print_genomes :
        for sp in args.print_genomes:
            entries = utils.get_specie_entries(blocks, sp)
            specie_genome = utils.thread_specie_genome(entries)
            utils.print_out_genome_thread(specie_genome, blocks)

