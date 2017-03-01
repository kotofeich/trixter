#!/usr/bin/env python

import argparse
from matplotlib import colors
from collections import defaultdict
from itertools import groupby
import matplotlib
matplotlib.use('Agg')
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

DIVISOR = 350000000.

#colorkey from etc/colors.ucsc.conf
colors = {'chr1': '#483d8b',
'chr2'  : '#87ceeb',
'chr3'  : '#3cb371',
'chr4'  : '#7fff00',
'chr5'  : '#ffff00',
'chr6'  : '#daa520',
'chr7'  : '#cd5c5c',
'chr8' : '#f4a460',
'chr9'  : '#b22222',
'chr10' : '#ffa07a',
'chr11' : '#ff8c00',
'chr12' : '#ff4500',
'chr13' : '#ff1493',
'chr14' : '#b03060',
'chr15' : '#9370db',
'chr16' : '#20b2aa',
'chr17' : '#000080',
'chr18' : '#778899',
'chr19' : '#ffe4c4',
'chr20' : '#556b2f',
'chr21' : '#4682b4',
'chr22' : '#4169e1',
'chr23' : '#7fff00',
'chrX' : '#cd5c5c'}

def parse_data(data,func=None):
    d = dict()
    current_chr=''
    with open(data) as f:
        for line in f:
            line = line.strip().split()
            if not line:
                current_chr = ''
            if len(line) == 1:
                if not current_chr:
                    d[line[0]] = []
                    current_chr = line[0]
                    first = True
                else:
                    if func:
                        line[0]=func(line[0])
                    d[current_chr].append([line[0]])
            elif len(line) == 3:
                val = map(int,line)
                if val[1]>val[2]:
                    start = val[2]
                    end = val[1]
                else:
                    start = val[1]
                    end = val[2]
                if first:
                    first = False
                    d[current_chr][-1].append(start)
                else:
                    d[current_chr][-1].append(end)
    return d
'''
def summarize(data):
    summed_data = defaultdict(list)
    for c in data:
        for k,g in groupby(data[c],lambda x: x[0]):
            s = sum(map(lambda x: x[1], g))/DIVISOR
            summed_data[c].append((k,s))
    return summed_data
   '''
    
def plot_everything(data,k):
    width = 0.03
    x_start = 0.02
    color = 'red'
    pp=PdfPages('figure1.pdf')
    fig = plt.figure(1, figsize=(9, 6))
    ax = fig.add_subplot(111)
    i = 0
    rectangles = []
    names = []
    for c in sorted(data.keys()):
        i += 1
        if i == 10:
            leg=plt.legend(rectangles,names,ncol=3, fontsize='xx-small', mode="expand", borderaxespad=0.)
            fig.patch.set_visible(False)
            ax.axis('off')
            plt.savefig(pp, format='pdf')
            plt.clf()
            plt.close()
            pp.close()
            rectangles = []
            names = []
            pp=PdfPages('figure2.pdf')
            fig = plt.figure(1, figsize=(9, 6))
            ax = fig.add_subplot(111)
            x_start = 0.02
        y_previous = 0
        r =  Rectangle((x_start+0.03, y_previous+0.04), width, k[c]/DIVISOR, facecolor = 'white')
        ax.add_patch(r)
        ax.annotate(c, (x_start+0.04, (y_previous+0.04+k[c]/DIVISOR)/2), color='black', weight='bold', 
                        fontsize=6, ha='center', va='center',rotation='vertical')
        x_start += 1.5*width 
        for block in data[c]:
            print block
            y_previous = block[1] / DIVISOR
            height = (block[2] - block[1]) / DIVISOR
            print block[0], colors[block[0]]
            r = Rectangle((x_start+0.03, y_previous + 0.04), width, height, facecolor = colors[block[0]],label=block[0],linewidth=0)
            if not block[0] in names:
                rectangles.append(r)
                names.append(block[0])
            ax.add_patch(r)
            y_previous += height
        x_start += 1.5*width 
    leg=plt.legend(rectangles,names,ncol=3, fontsize='xx-small', mode="expand", borderaxespad=0.)
    fig.patch.set_visible(False)
    ax.axis('off')
    plt.savefig(pp, format='pdf')
    pp.close()
    
 
def parse_karyotype(k):
    kar = {}
    with open(k) as f:
        for line in f:
            line = line.strip().split()
            kar[line[0]] = int(line[1])
    return kar

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='genome1000K')
    parser.add_argument('k', help='reference karyotype obtained with faSize')
    args = parser.parse_args()

    data = parse_data(args.file, lambda x: x.split('HomoSapiens.')[1])
    k = parse_karyotype(args.k)
    for d in data:
        print d, data[d]
    #data = summarize(data)
    plot_everything(data,k)

