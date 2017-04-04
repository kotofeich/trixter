#!/usr/bin/env python

import argparse
from time import gmtime, strftime
import model
import utils
import rearrangements_type

def log_time(log):
    log.write(strftime("%Y-%m-%d %H:%M:%S", gmtime())+'\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('maf', help='maf file')
    parser.add_argument('--species', nargs='+')
    args = parser.parse_args()
    log = open('./maf_rear.scaffold7.log','w+')
    if not args.species or len(args.species) != 2:
        exit('Number of species must be 2')
    log_time(log)
    log.write('parsing maf...\n')
    blocks = model.parse_maf(args.maf)
    log_time(log)
    log.write('collecting species entries...\n')
    entries = utils.get_specie_entries(blocks, args.species[0])
    log_time(log)
    log.write('threading specie genome...\n')
    ref_genome = utils.thread_specie_genome(entries)
    #print ref_genome
    log_time(log)
    log.write('collecting species entries...\n')
    entries = utils.get_specie_entries(blocks, args.species[1])
    log_time(log)
    log.write('threading specie genome...\n')
    target_genome = utils.thread_specie_genome(entries)
    #print target_genome
    entries = []
    log_time(log)
    log.write('grouping by reference..\n')
    target_genome_grouped = utils.group_by_ref(ref_genome, target_genome)
    '''
    for e in target_genome_grouped:
        for k in e:
            k.print_out()
    print
    for e in ref_genome:
        for k in e:
            k.print_out()
    '''
    log_time(log)
    log.write('printing genome thread..\n')
    path='./FCscaffold7.'+args.species[0]
    utils.print_out_genome_thread(ref_genome,path)
    log_time(log)
    log.write('printing genome thread..\n')
    path='./FCscaffold7.'+args.species[1]
    utils.print_out_genome_thread(target_genome,path)
    ##ref_genome, target_genome_grouped = utils.normalize(ref_genome, target_genome_grouped)
    log_time(log)
    log.write('analyzing rearrangements..\n')
    for c in target_genome_grouped:
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
            #start is the block before the trasposition
            if not this_prev in map(lambda x: x[1], trp):
                to_start = this_prev
            if not this_next in map(lambda x: x[1], trp):
                to_end = this_next
            if to_start != -1 and to_end != -2: 
                #from_start = utils.find_prev_block_in_specie(this_trp[0],specie2)
                #from_end = utils.find_next_block_in_specie(this_trp[-1], specie2)
                #l = str(from_start.block_id)+'\t' if from_start else 'None\t'
                #l += str(from_end.block_id)+'\t' if from_end else 'None\t'
                #l += str(to_start.block_id)+'\t' if to_start else 'None\t'
                #l += str(to_end.block_id) if to_end else 'None'
                #l='transposition'
                print 'transposition'
                for t in this_trp:
                    #l += '\t'+str(t.block_id)
                    t.print_out_short()
                #print l
                to_start = -1
                to_end = -2
                this_trp = []
    
        main_chrom, trl = rearrangements_type.check_translocations(c)
        for e in trl:
            print 'translocation: from chromosome', main_chrom
            for x in e:
                x.print_out_short()
        #if trl:
        #    print 'overall translocations:', len(trl)
