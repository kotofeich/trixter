#!/usr/bin/env python

import sys

if __name__ == '__main__':
    file_name = sys.argv[1]
    summary = []
    with open(file_name) as f:
        name1 = None
        name2 = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if 'translocation:' in line:
                if name1 and name2:
                    summary.append(sorted([name1, name2]))
                name1 = line.split("from chromosome ")[1]
            elif 'seq_id: ' in line:
                name2 = line.split('seq_id: ')[1]
                name2 = name2.split(' block_id: ')[0]
                name2 = name2.split('.')[1]
    summary = sorted(summary, key=lambda x:x[0])
    for x in summary:
        print x[0], '+', x[1]
