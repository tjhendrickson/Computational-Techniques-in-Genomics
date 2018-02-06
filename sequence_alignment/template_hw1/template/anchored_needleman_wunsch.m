function [score, alignment1, aligner, alignment2] = anchored_needleman_wunsch(seq1, seq2, match, mismatch, gap, matched_regions)
% anchored_needleman_wunsch Perform the anchored version of the 
% Needleman-Wunsch algorithm.
%
% Input variables:
% seq1: string containing the first sequence
% seq2: string containing the second sequence
% match: numeric value which represents the score for match
% mismatch: numeric value which represents the penalty for mismatch
% gap: numeric value which represents the penalty for gap
% matched_regions: contains regions which are know to be aligned. The
% regions may contain matches and/or mismatches. The first column contains
% the starting position of the first sequence; the second column contains
% the ending position of the first sequence; the third column contains the
% starting position of the second sequence; and the forth column contains
% the ending position of the second sequence. E.g.: a row containing [17,
% 20, 76, 79] means that seq1(17:20) should be aligned to seq2(76:79),
% containing only matches and/or mismatches in this region.
%
% Output variables:
% score: numeric value containing the alignment score
% alignment1: string containing the alignment (including gaps) of the first
% sequence
% aligner: string containing symbols (pipe |) which aligns the 2 sequences
% alignment2: string containing the alignment (including gaps) of the
% second sequence
%
% Output example:
% AV--TNAGQLV---S  <--- alignment1
% ||  || |||    |  <--- aligner
% AQVSTN-GQL-AQVT  <--- alignment2
%
%
%%%%%%%%%%%%%% YOUR CODE STARTS HERE

%first create all sub problems (i.e. sequences prior to matched regions,
%sequences in matched regions, and sequences after matched regions,
%and place it within a matrix
sub_problems=[ 1 matched_regions(1,1)-1, 1, matched_regions(1,3)-1;
matched_regions(1,:);
matched_regions(1,2)+1 matched_regions(2,1)-1 matched_regions(1,4)+1 matched_regions(2,3)-1;
matched_regions(2,:);
matched_regions(2,2)+1 matched_regions(3,1)-1 matched_regions(2,4)+1 matched_regions(3,3)-1;
matched_regions(3,:);
matched_regions(3,2)+1 length(seq1) matched_regions(3,4)+1 length(seq2)
];

%initialize score, alignment, aligner, and alignment2
score=0;
%define alignment1, aligner, alignment2
alignment1=['']; %corresponds to sequence 1
aligner=[''];
alignment2=['']; %corresponds to sequence 2 

%We are interested in iterating through the number of rows
for i=1:length(sub_problems)
    
    %subsequence one corresponds to the ith row and column 1 for start and
    %ith row and column 2 to end
    subseq1=seq1(sub_problems(i,1):sub_problems(i,2));
    
    %ith row column 3 = start, ith row column 4 = end
    subseq2=seq2(sub_problems(i,3):sub_problems(i,4));
    
     
    %call needleman wunsch
    
    [score_need_wunsch,alignment1_need_wunsch,aligner_need_wunsch,alignment2_need_wunsch] = needleman_wunsch(subseq1, subseq2, match, mismatch, gap,matched_regions);
    
    %report back score, alignment1, aligner, and alignment2
    score=score+score_need_wunsch; 
    alignment1=[alignment1;alignment1_need_wunsch];
    aligner=[aligner;aligner_need_wunsch]; 
    alignment2=[alignment2;alignment2_need_wunsch];

end


end

