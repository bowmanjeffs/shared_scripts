# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 10:50:36 2013

@author: Jeff

Takes as input a single aligned fasta file.  be sure to set appropriate gap
characters and the minimum number of positions.

The script will remove all positions before the last start position and after
the earliest end position.  it will 

python filter_seqs_2.py <file in>

v5 allows sequences in alignment to take up multiple lines

"""
import re
import sys

file_in = sys.argv[1]
#file_in = 'DUF900_pro_aligned.fasta'  ## use this for testing
name = re.sub('.fasta', '', file_in)
log = open(name+'.filter.log', 'w')

gap_character = ['-', '.']
min_positions = 500 ## minimum positions allowed after filtering

for character in gap_character:
    print >> log, 'gap = '+character
    
print >> log, '\n',

## find start and end for each sequence
                    
start_end = {}
max_starts = []
min_ends = []
starts = set()
ends = set()

print 'filtering to last start and first end'
print >> log, 'filtering to last start and first end'

l = 0

with open(file_in, 'r') as fasta_in:
    lines = ''
    for line in fasta_in:
        
        if line.startswith('>') == False:
            line = line.rstrip()
            lines = lines + line
            
        else:
            l = l + 1
            if l != 1:
                rlines = lines[::-1]
                
                l = l + 1
                
                for i,p in enumerate(lines):
                    if p not in gap_character:
                        start = i - 1
                        break
                        
                for i,p in enumerate(rlines):
                    if p not in gap_character:
                        end = len(rlines) - i
                        break
                                
                print l, start, end
                print >> log, l, start, end
                start_end[l] = start, end
                starts.add(start)
                ends.add(end)
                lines = ''

## make sure you get the last line!

l = l + 1

for i,p in enumerate(lines):
    if p not in gap_character:
        start = i - 1
        break
        
for i,p in enumerate(rlines):
    if p not in gap_character:
        end = len(rlines) - i
        break
                
print l, start, end
print >> log, l, start, end
start_end[l] = start, end
starts.add(start)
ends.add(end)
            
min_end = min(ends)
max_start = max(starts)   
print >> log, 'max start / min end =', max_start, '/', min_end 

## eliminate late starts and early finishes until minimum alignment length reached

bad = []

n = 0
starts = sorted(starts, reverse = True)
ends = sorted(ends)

while(min_end - max_start) < min_positions:
    print 'min_end - max_start =', min_end - max_start
    print >> log, 'max start - min end is below cutoff, eliminating some sequences'        
    for key in start_end.keys():
        if key not in bad:
            if start_end[key][0] == max_start:
                bad.append(key)
            elif start_end[key][1] == min_end:
                bad.append(key)
                
    min_end = ends[n]
    max_start = starts[n]
    n = n + 1
 
## look for gap only columns among good sequences

gap = {}
       
with open(file_in, 'r') as fasta_in:
    nseq = 0
    lines = ''
    for line in fasta_in:
        
        if line.startswith('>') == False:
            line = line.rstrip()
            lines = lines + line
        
        else:
            nseq = nseq + 1
            if nseq != 1:
            
                nseq = nseq + 1
                if nseq not in bad:
                    print nseq
                    
                    for i,p in enumerate(lines):
                        if p in gap_character:
                            try:
                                temp = gap[i]
                                temp = temp + 1
                                gap[i] = temp
                            except KeyError:
                                gap[i] = 1
                lines = ''
                
nseq = nseq + 1
if nseq != 1:

    nseq = nseq + 1
    if nseq not in bad:
        print nseq
        
        for i,p in enumerate(lines):
            if p in gap_character:
                try:
                    temp = gap[i]
                    temp = temp + 1
                    gap[i] = temp
                except KeyError:
                    gap[i] = 1
                            
## generate filter and remove gaps
                    
flter = {}

for key in range(0, len(lines)):
    if key < max_start:
        flter[key] = 0
    elif key > min_end:
        flter[key] = 0
    elif key in gap.keys():
        if gap[key] == nseq - len(bad): 
            flter[key] = 0
        else:
            flter[key] = 1
    else:
        flter[key] = 1

flter_list = []        
for key in flter.keys():
    if flter[key] == 1:
        flter_list.append(key)
            
with open(file_in, 'r') as fasta_in, open(name+'.filter.fasta', 'w') as fasta_out:
    l = 0
    seq = ''
    for line in fasta_in:
        if line.startswith('>'):
            if l != 0:
                seq_out = ''
                for i,p in enumerate(seq):
                    if flter[i] == 1:
                        seq_out = seq_out + p
                print >> fasta_out, seq_out
            l = l + 1
            seq = ''
            if l not in bad:
                print >> fasta_out, line,
            else:
                line = line.strip()
                print >> log, line+' is bad'
        
        else:
            if l not in bad:
                line = line.rstrip('\n')
                seq = seq + line
    print >> fasta_out, seq_out
                
print >> log, '\n'
print >> log, 'filter:'
            
for key in flter.keys():
    print >> log, flter[key],

log.close()
                