Scope of this document - To provide a description on how to run the needleman wunsch algorithm and anchored needleman wunsch algorithm in Matlab

-The first step is to download MATLAB version 2016a. The user is recommended to download either this version or a more recent version. If earlier than 2016a it is possible that syntax or functions may have changed or did not exist. 

-Once MATLAB has been downloaded the user can then start working on setting up directory/folder intrastructure to run the alignments. If not already downloaded, download the matlab scripts anchored_needleman_wunsch.m, needleman_wunsch.m, and main.m along with the fasta files that you wish to sequence, and matched regions to be used with the anchored needleman wunsch. Once you have all of this material place all of the files into the same folder, the scripts will not work if the files are in different folders/directories. 

-In order to set up the scripts to run all the user has to do is open main.m and modify a few variables at the beginning of the file. The variables to be modified are: seq1, seq2, proteinAligned, and matched_regions. This is pretty self explanatory but the two sequences to be aligned go within each of the seq variables. Next, place a string of the gene or protein that one is aligning. Finally, place the string corresponding to the matched regions text file within th matched_regions variable.  

-Now the script is ready to run! The script should run through without problem and output a variety of data. First the alignment score, and alignment for both the needleman wunsch and anchored needleman wunsch algorithms will be appended to a text file named sequence_analysis.txt. Later a histogram will be generated with the output of 10,000 random sequence alignments along the actual sequence.

For questions or comments please contact Tim Hendrickson at hendr522@umn.edu.  


