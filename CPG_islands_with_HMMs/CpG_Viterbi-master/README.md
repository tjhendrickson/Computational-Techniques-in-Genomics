CpGviterbi
==========
CpGviterbi algorithm which creates a hidden markov model to detect CpG islands within an input DNA Sequence. It uses the viterbi algorithm to find the most probable state path among the sequence.

You can call the function in matlab by:

*CpGviterbi('seq1.fasta')*

This will return the found CpG islands, along with the number of CpG islands, the number of down-stream genes that begin coding within 500bps of a CpG island, and a matrix of the lengths of the CpG islands.

Along with these outputs, each CpG island is printed out in matlab in the provided format:

- CpG island 1: 348 bp (739 - 1086) gene_name ; number_of_motif_hits
- CpG island 2: 1148 bp (1508 - 2655) ; number_of_motif_hits
- CpG island 3: 971 bp (4029 - 4999) gene_name ; number_of_motif_hits
- Total CpG islands found: 3 ; 2 out of 3 islands are followed by a coding region 

Note: The CpG islands are displayed in reverse order because of the traceback algorithm (e.g. CpG Island 1 is the last CpG island in the sequence).

The motif detection loop is commented out for speed. Uncomment it to display them (line 167).


Other files in directory:
- motif_hits.txt -- map of the motifs each with their corresponding CpG island
- results.txt -- The CpG islands listed with the format above
- CpGviterbi.m -- source code
