#!/usr/bin/env python
import gzip
import numpy as np
import bioinfo
import argparse
from itertools import permutations

def get_args():
    parser = argparse.ArgumentParser(description="Demultiplex")
    parser.add_argument("-r1", "--read1", help="R1 File name", required=True, type=str)
    parser.add_argument("-r2", "--read2", help="R2 File name", required=True, type=str)
    parser.add_argument("-r3", "--read3", help="R3 File name", required=True, type=str)
    parser.add_argument("-r4", "--read4", help="R4 File name", required=True, type=str)
    parser.add_argument("-i", "--index", help="Known Indexes File name", required=True, type=str)
    return parser.parse_args()
args = get_args()

read1 = gzip.open(args.read1, "rt")
index1 = gzip.open(args.read2, "rt")
index2 = gzip.open(args.read3, "rt")
read2 = gzip.open(args.read4, "rt")
unknown_r1 = gzip.open("Unknown_R1.fastq.gz", 'wt')
unknown_r2 = gzip.open("Unknown_R2.fastq.gz", 'wt')
unmatched_r1 = gzip.open("Unmatched_R1.fastq.gz", 'wt')
unmatched_r2 = gzip.open("Unmatched_R2.fastq.gz", 'wt')
#Create set of 24 known indexes (forward/given)
indexes = set()
#rev_indexes = set()
with open(args.index, "r") as indexes_file:
    for i in indexes_file:
        i = i.strip()
        indexes.add(i)
#         rev = bioinfo.rev_comp(i)
#         rev_indexes.add(rev)
# # print(indexes)
# print(rev_indexes)

# Create counters to keep track of number of reads going into each file: matched, unmatched, unknown
unmatched_counter = 0
unknown_counter = 0
low_qual = 0
matched_counter = 0
barcode_matched= {}
# Create empty dictionary to keep track of the amount of times each barcode had went to the matched file: barcodes_matched
matched_r1={}
matched_r2={}
for i in indexes:
    matched_r1[i] = gzip.open(i+"_R1.fastq.gz", "wt") #w or a+
    matched_r2[i] = gzip.open(i+"_R2.fastq.gz", "wt")

#-create dictionary of possible permutations
perms = permutations(indexes, 2)
known_hopped = {}
for i, j in perms:
    key = i+'_'+j
    known_hopped[key] = 0
#print(known_matched.items())
#-if(index1_revcomp2 in dictionary) --> increment
#Open all 4 files and read in parallel (will be looping through the records)
#while true:
    #loop through every record
    #Collect all 4 lines of the record for each file and place into 4 temporary arrays of size 4 (will save each line as a position in the array)

counter = 0 #Keep a counter called to keep track of the total number of lines in the file.
while True:
    r1_line1 = read1.readline()
    i1_line1 = index1.readline()
    i2_line1 = index2.readline()
    r2_line1 = read2.readline()
    counter+=1
    if(r1_line1 == ""): #all files should have the same number of records, so only need to check one of the lines
        counter-=1 #reached the end of the file
        break #reach the end of the file, break out of the while loop
    elif(r1_line1.startswith('@')):
        read1_rec = np.array([r1_line1.strip(), read1.readline().strip(), read1.readline().strip(), read1.readline().strip()])
        index1_rec = np.array([i1_line1.strip(), index1.readline().strip(), index1.readline().strip(), index1.readline().strip()])
        index2_rec = np.array([i2_line1.strip(), index2.readline().strip(), index2.readline().strip(), index2.readline().strip()])
        read2_rec = np.array([r2_line1.strip(), read2.readline().strip(), read2.readline().strip(), read2.readline().strip()])       
        counter+=3
        if(counter/4 % 100000 == 0):
            print(counter/4, "record")

        #START DEMULTIPLEXING HERE:

        #Read the sequence line of the nth record for all 4 files
        #For each record in the arrays:
        revcomp_ind2 = bioinfo.rev_comp(index2_rec[1]) #save the reverse comp of index2 to a variable
        # rev_ind2 = ""
        # for i in reversed(index2_rec[1]):
        #     rev_ind2+=i
        #Check if any N's in the index 1 or index 2 sequences. If either containts an N --> unknown file, increment unknown counter
        if ('N' in index1_rec[1]) or ('N' in index2_rec[1]):
            bioinfo.write_output(unknown_r1, read1_rec, index1_rec[1], revcomp_ind2)
            bioinfo.write_output(unknown_r2, read2_rec, index1_rec[1], revcomp_ind2)
            unknown_counter +=1
            #print('N')
            continue #go back to the top of the while loop and start checking the next record
        #if any base pair position has Qscore < 30 for either index --> unknown file, increment unknown counter
        for ind, value in enumerate(index1_rec[3]):
            qscore_1 = bioinfo.convert_phred(value)
            qscore_2 = bioinfo.convert_phred(index2_rec[3][ind])
            bool_break = False
            if(qscore_1 < 30) or (qscore_2 < 30):
                bioinfo.write_output(unknown_r1, read1_rec, index1_rec[1], revcomp_ind2)
                bioinfo.write_output(unknown_r2, read2_rec, index1_rec[1], revcomp_ind2)
                unknown_counter +=1
                low_qual+=1
             #   print('Qscore')
                bool_break = True
                break #break out of the for loop, don't need to keep checking the rest of the Qscores because one was less than 30
        if(bool_break == True):
            continue #go back to the top of the while loop and start checking the next record
        #Check if index 1 in set:
        #if index 1 is not in the set --> unknown file, increment unknown counter
        # print(indexes)
        # print(rev_indexes)
        if(index1_rec[1] not in indexes) or (revcomp_ind2 not in indexes):
            bioinfo.write_output(unknown_r1, read1_rec, index1_rec[1], revcomp_ind2)
            bioinfo.write_output(unknown_r2, read2_rec, index1_rec[1], revcomp_ind2)
            unknown_counter +=1
            #print('unknown index', index1_rec[1], revcomp_ind2, sep = " ")
            continue ##go back to the top of the while loop and start checking the next record
        #if both indexes in set --> continue filtering (check if they match)
        else:
        #elif(index1_rec[1] in indexes) and (revcomp_ind2 in rev_indexes):
            #Check if matching: does index 1 == rev comp index 2
            #No --> unmatched file, increment unmatched counter
            
            if(index1_rec[1] != revcomp_ind2):
                bioinfo.write_output(unmatched_r1, read1_rec, index1_rec[1], revcomp_ind2)
                bioinfo.write_output(unmatched_r2, read2_rec, index1_rec[1], revcomp_ind2)
                unmatched_counter +=1
                unmatched_key = index1_rec[1] + '_' + revcomp_ind2
                known_hopped[unmatched_key] +=1
                #print('unmatched', index1_rec[1], revcomp_ind2)
                continue #go back to the top of the while loop and start checking the next record
            #Yes --> index1_rcindex2_Read#.fq (matched files) (rcindex2 = reverse complement of ind2), increment matched counter
            #also will have barcodes_matched[key] = value, 
                #if key is in the dictionary, increment the value
                #if key is not in the dictionary, add it and set equal to 1
            elif(index1_rec[1] == revcomp_ind2):
                #print("matched", index1_rec[1], revcomp_ind2)
                matched_counter+=1
                if(index1_rec[1] not in barcode_matched):
                    barcode_matched[index1_rec[1]] = 1
                    bioinfo.write_output(matched_r1[index1_rec[1]], read1_rec, index1_rec[1], revcomp_ind2)
                    bioinfo.write_output(matched_r2[index1_rec[1]], read2_rec, index1_rec[1], revcomp_ind2)
                    #print('added to dict')
                else:
                    #print('increment')
                    barcode_matched[index1_rec[1]] += 1
                    #print(barcode_matched)
                    bioinfo.write_output(matched_r1[index1_rec[1]], read1_rec, index1_rec[1], revcomp_ind2)
                    bioinfo.write_output(matched_r2[index1_rec[1]], read2_rec, index1_rec[1], revcomp_ind2)

for i in indexes:
    matched_r1[i].close()
    matched_r2[i].close()
read1.close()
index1.close()
index2.close()
read2.close()
unknown_r1.close()
unknown_r2.close()
unmatched_r1.close()
unmatched_r2.close()

print('Number of Records in unknown: '+str(unknown_counter))
print('Number of low quality records: '+str(low_qual))
print('Number of unmatched: '+str(unmatched_counter))
print('Hopped Dictionary:')
for key,value in known_hopped.items():
    print(key,value)
print('Number of matched: '+str(matched_counter))
print('Matched dictionary:')
for key,value in barcode_matched.items():
    print(key,value)

                        