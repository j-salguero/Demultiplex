#!/usr/bin/env python

import gzip
import numpy as np
import bioinfo
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Demultiplex")
    parser.add_argument("-f", "--file", help="File name", required=True, type=str)
    parser.add_argument("-l", "--length", help="Seq length", required=True, type=int)
    parser.add_argument("-g", "--graph", help="Graph title", required=True, type=str)
    return parser.parse_args()
args = get_args()

num_lines = 0
#loop through every record
with gzip.open(args.file, "rt") as opened: #rt = read text mode
    #Make the empty list
    #each position is a pos, with an empty list where I'll put the qscores at that position
    all_qscores: list = []
    all_qscores = bioinfo.init_list(all_qscores, args.length)
    i = 0 #line number counter
    for line in opened:
        i+=1
        num_lines+=1
        if i%4 == 0: #if it's a sequence line
            #line = line.decode("ascii") #turn the zipped bits into the ascii characters
            line = line.strip()
            rec_num= i//4 #record number
            if rec_num%100000 == 0:
                print(f'curr record {rec_num}')
    #convert the Phred quality score from a letter to its corresponding number 
    #and add it to an ongoing sum inside your list of the quality scores for each base pair.
            #counter_phred = 0
            for ind, value in enumerate(line):
                phred_pos = bioinfo.convert_phred(value)
                all_qscores[ind] +=phred_pos
                # all_qscores[ind].append(phred_pos)

mean_qscores = np.zeros(args.length,dtype=np.float64)
for ind, value in enumerate(all_qscores):
    mean_qscores[ind] = value/(num_lines//4)

print(mean_qscores)
import matplotlib.pyplot as plt
#%matplotlib inline
plt.bar(range(args.length), mean_qscores)
plt.xlabel("Base Position")
plt.ylabel("Mean Qscore")
plt.title(args.graph)
plt.savefig(args.graph+'.png')