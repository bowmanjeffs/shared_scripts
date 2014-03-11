# -*- coding: utf-8 -*-
"""
Created on Tue Mar 11 09:16:05 2014

@author: Jeff
"""
p = 30 # phred score cutoff
m = 15 # max bases to trim off front and back

import sys

ffasta = sys.argv[1]
rfasta = sys.argv[2]
ifasta = sys.argv[3]

#ffasta = 'small.fastq'
#rfasta = 'rsmall.fastq'
#ifasta = 'ismall.fastq'

bad = set()

from Bio import SeqIO

with open('log.txt', 'w') as log:
    n = 0
    for record in SeqIO.parse(ffasta, 'fastq'):
        i = 0
        n = n + 1
        print n
        init_l = len(record.letter_annotations["phred_quality"])
        init_mean = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])
        while sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"]) < p:
            i = i + 1
            if i == m:
                bad.add(record.id)
                break
            else:
                record = record[1:-1]
        print >> log, 'forward', record.id, init_l, len(record.letter_annotations["phred_quality"]), init_mean, sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])
    
    n = 0        
    for record in SeqIO.parse(rfasta, 'fastq'):
        i = 0
        n = n + 1
        print 'reverse', n
        if record.id not in bad:
            init_l = len(record.letter_annotations["phred_quality"])
            init_mean = sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])
            while sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"]) < p:
                i = i + 1
                if i == m:
                    bad.add(record.id)
                    break
                else:
                    record = record[1:-1]
            print >> log, 'reverse', record.id, init_l, len(record.letter_annotations["phred_quality"]), init_mean, sum(record.letter_annotations["phred_quality"]) / len(record.letter_annotations["phred_quality"])
            
print len(bad), 'discarded out of', n

print 'writing forward fasta'
with open(ffasta+'.trim', 'w') as fasta_out:
    for record in SeqIO.parse(ffasta, 'fastq'):
        if record.id not in bad:
            SeqIO.write(record, fasta_out, 'fastq')

print 'writing reverse fasta'            
with open(rfasta+'.trim', 'w') as rfasta_out:
    for record in SeqIO.parse(rfasta, 'fastq'):
        if record.id not in bad:
            SeqIO.write(record, rfasta_out, 'fastq')

print 'writing index fasta'            
with open(ifasta+'.trim', 'w') as index_out:
    for record in SeqIO.parse(ifasta, 'fastq'):
        if record.id not in bad:
            SeqIO.write(record, index_out, 'fastq')
