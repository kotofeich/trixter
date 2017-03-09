#!/usr/bin/env python

import argparse
import model

def print_out(blocks, ref, target):
    for b in blocks:
        ref_data = ''
        target_data = ''
        dup_alarm = False
        for e in b.entries:
            if e.get_specie() == ref:
                if ref_data:
                    dup_alarm =True
                    break
                ref_data = '\t'.join([e.get_chrom(),str(e.start),str(e.end)]) 
            if e.get_specie() == target:
                if target_data:
                    dup_alarm =True
                    break
                marker_name = target+'.'+str(e.block_id)
                target_data = '\t'.join([e.get_chrom(),str(e.start),str(e.end),marker_name,marker_name,marker_name]) 
        if dup_alarm:
            dup_alarm = False
            continue
        elif ref_data and target_data:
            print '\t'.join([target_data,ref_data])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='blocks_coords.txt')
    parser.add_argument('ref', help='reference genome name')
    parser.add_argument('target', help='genome name')
    args = parser.parse_args()
    blocks = model.parse_blocks(args.file)
    print_out(blocks,args.ref,args.target)
