function [score, alignment1, aligner, alignment2] = needleman_wunsch(seq1, seq2, match, mismatch, gap,matched_regions)
% needleman_wunsch Perform the Needleman-Wunsch algorithm.
%
% Input variables:
% seq1: string contalignment1(i)ning the first sequence
% seq2: string contalignment1(i)ning the second sequence
% match: numeric value which represents the score for match
% mismatch: numeric value which represents the penalty for mismatch
% gap: numeric value which represents the penalty for gap
%
% Output variables:
% score: numeric value contalignment1(i)ning the alignment score
% alignment1: string contalignment1(i)ning the alignment (including gaps) of the first
% sequence
% aligner: string contalignment1(i)ning symbols (pipe |) which aligns the 2 sequences
% alignment2: string contalignment1(i)ning the alignment (including gaps) of the
% second sequence
%
% Output example:
% AV--TNAGQLV---S  <--- alignment1
% ||  || |||    |  <--- aligner
% AQVSTN-GQL-AQVT  <--- alignment2
%
%
%%%%%%%%%%%%%% YOUR CODE STARTS HERE

%create a MxN matrix with M corresponding to seq2 and N corresponding to seq1

%seq1 = number of columns

%seq2 = number of rows

%scoring matrix
F=zeros(length(seq1)+2,length(seq2)+2); %which results in a matrix 232x382

%traceback matrix
Ptr=F;

%now define number of rows as rowsF and columns as columns F
[rowsF, columnsF]=size(F);

%set position at intersection between first row and column of matrix to be zero
F(2,2)=0;

%set up numbers corresponding to seq1 and seq2 locations
F(3:end,1)=(1:length(seq1));
F(1,3:end)=(1:length(seq2));

%add gap scores starting with second rows and second column to matrix 
 F(3:end,2)=(1:length(seq1))*gap;
 F(2,3:end)=(1:length(seq2))*gap;

%set up the double nested for loop to iterate through the two sequences to
%be compared
for i=3:rowsF; % iterate from 1 to length of seq1
    for j=3:columnsF; %iterate from 1 to length of seq2
        
        %first determine if there is a match between seq1(i) and seq2(j)
        if strcmp(seq1(F(i,1)),seq2(F(1,j)))
            s=match;
        else
            s=mismatch;
        end
        
        %now score the three possibilities
        
        %if previous location was a match or mismatch 
        diagonal=F(i-1,j-1)+s;
        
        %if previous location was a gap from left
        left=F(i,j-1)+gap;
        
        %if previous location was a gap from up
        up=F(i-1,j)+gap;
        
        %put all three scores into a row vector
        outcomes=[left,diagonal,up];
        
        %now determine max of outcomes, which is the optimal sub alignment
        F(i,j)=max(outcomes);
        
        % what is the traceback for this optimal sub alignment with the key
        % 1 = left, 2 = diagonal, 3 = up
        if F(i,j)==outcomes(1);
            Ptr(i,j)=1;
        elseif F(i,j)==outcomes(2);
            Ptr(i,j)=2;
        else
            Ptr(i,j)=3;
        end
        
    
    end
end

score=F(rowsF,columnsF);
%define alignment1, aligner, alignment2
alignment1=['']; %corresponds to sequence 1
aligner=[''];
alignment2=['']; %corresponds to sequence 2 

%now take Ptr and traceback maximum score

%how many rows in matrix
 i=rowsF;
%how many columns in matrix
 j=columnsF;
 %while sequence location is above 1 in seq1 or seq2 keep running otherwise
 %quit
 while F(i,1)>1 || F(1,j)>1
  %while i>=2 || j>   
        %if position Ptr(i,j) in traceback = 1 place a gap in alignment1
        %fill in the amino acid in alignment 2
        %place a space for the aligner
        %then subtract j by one
         if Ptr(i,j)==1 && F(1,j)>=1;
             alignment1=[seq1(F(i,1));alignment1];
             alignment2=['-';alignment2];
             aligner=[' ';aligner];
             i=i-1;
                      
         %if position Ptr(i,j) in trackback = 2 retrieve the amino acid for
         %alignment 1 and 2 and place a '|' in aligner
         %subtract i and j by one
         elseif Ptr(i,j)==2 && F(i,1)>=1 && F(1,j)>=1
             alignment1=[seq1(F(i,1));alignment1];
             alignment2=[seq2(F(1,j));alignment2];
             aligner=['|';aligner];
             j=j-1;
             i=i-1;
             
          %if position Ptr(i,j) in traceback = 3 place a gap in alignment2
        %fill in the amino acid in alignment 1
        %place a space for the aligner
        %then subtract i by one    
         else
             alignment1=['-';alignment1];
             alignment2=[seq2(F(1,j));alignment2];
             aligner=[' ';aligner];
             j=j-1;
         end
        
         %create a catch in case j runs to 0 in order to place gaps in
         %second alignmet and AA's in first alignment
         
         if F(1,j)==0
             while F(i,1)>=1
                 alignment1=[seq1(F(i,1));alignment1];
                 alignment2=['-';alignment2];
                 aligner=[' ';aligner];
                 i=i-1;
             end
         end
 
end
