shared_scripts
==============
filter_seqs_4.py
--------------
This script takes an aligned fasta file and returns a fast file with gap only columns removed, and all positions preceding the last start and following the earliest finish removed.  If the resulting alignment is below a specified length sequences are removed from the alignment and the filter recalculated until a satisfactory length is reached.