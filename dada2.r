#!/usr/bin/Rscript

# See tutorial at https://benjjneb.github.io/dada2/tutorial.html.  This script will
# QC and identify valid unique reads, then assemble.  It will execute on all files in
# the directory "multiplexed", and create the directories "filtered", and "merged".
# The Python script deunique_dada2.py should be used to inflate all of the output tables
# into redundant fasta files for paprica.

library(dada2)

path <- 'demultiplexed'

fnFs <- sort(list.files(path, pattern = '-R1.fastq', full.names = T))
fnRs <- sort(list.files(path, pattern = '-R2.fastq', full.names = T))

sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)

pdf('quality_profiles.pdf', width = 6, height = 6)

for(i in 1:length(fnFs)){
	plotQualityProfile(fnFs[i])
	plotQualityProfile(fnRs[i])
	}
	
dev.off()

file_path <- file.path("filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(file_path, paste0(sample.names, "_R1.filt.fastq"))
filtRs <- file.path(file_path, paste0(sample.names, "_R2.filt.fastq"))

## multithreading only useful if multiple fastq files

out <- filterAndTrim(fnFs,
	filtFs,
	fnRs,
	filtRs,
	#multithread = T,
	minQ = 20,
	verbose = T)
	
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

pdf('error_rates.pdf', width = 6, height = 6)
plotErrors(errF, nominalQ = T)
plotErrors(errR, nominalQ = T)
dev.off()

derepFs <- derepFastq(filtFs, verbose = T)
derepRs <- derepFastq(filtRs, verbose = T)

## wouldn't it make more sense to assemble here?

names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errF, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

dir.create('merged')

for(name in names(mergers)){
	write.csv(mergers[[name]], paste0('merged/', name, '.csv'), quote = F, row.names = F)
	}
