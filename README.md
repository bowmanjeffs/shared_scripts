shared_scripts
==============
filter_seqs_4.py
--------------
This script takes an aligned fasta file and returns a fast file with gap only columns removed, and all positions preceding the last start and following the earliest finish removed.  If the resulting alignment is below a specified length sequences are removed from the alignment and the filter recalculated until a satisfactory length is reached.

python filter_seqs_4.py [aligned fasta]

trim_fastq_v2.py
--------------
This script is a lightweight trimmer.  It takes in a forward fastq, reverse fastq, and index fastq from a barcoded Illumina run and trims to a specified mean quality score.  Requires biopython.

python trim_fastq_v2.py [forward fastq] [reverse fastq] [index fastq]