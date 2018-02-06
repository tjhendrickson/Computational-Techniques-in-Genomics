%function [cpgIslands, cpgCount, codingCount, cpgLengthMtr] = CpGviterbi(seq)
%implementation of Hidden Markov Model for CpG Islands
%which uses the viterbi algorithm to find the best path

%toy samples
%[head, toy1] = fastaread('Toy_seq\seq1.fa');
%[head, toy2] = fastaread('Toy_seq\seq2.fa');
%input sequence
FASTAfilename='seq1.out';
%[head, seqx] = fastaread(seq);

fidIn = fopen(FASTAfilename, 'r');
seqx = textscan(fidIn, '%s');
seqx = strjoin(seqx{1}','');
fclose(fidIn);

%hmm generation values
p = 0.98;
q = 0.999;
%values for toy sequences
%p = 0.9999;
%q = 0.95;

%transition matrix for CpG islands 
%       |  A+     |   C+    |     G+      |    T+     |     A-     |    C-      |   G-     |    T-   |
TRANS = [0.180*p,   0.274*p,    0.426*p,    0.120*p,    ((1-p)/4),  ((1-p)/4),  ((1-p)/4),  ((1-p)/4) ;
         0.171*p,	0.368*p,	0.274*p, 	0.188*p,  	((1-p)/4),	((1-p)/4),	((1-p)/4),	((1-p)/4) ;
    	 0.161*p,	0.339*p,   	0.375*p,	0.125*p,    ((1-p)/4),	((1-p)/4),	((1-p)/4),	((1-p)/4) ;
    	 0.079*p,	0.355*p, 	0.384*p, 	0.182*p, 	((1-p)/4),	((1-p)/4),	((1-p)/4),	((1-p)/4) ;
    	((1-q)/4),	((1-q)/4),	((1-q)/4),	((1-q)/4),	0.300*q,  	0.205*q, 	0.285*q, 	0.210*q   ;
    	((1-q)/4),	((1-q)/4),	((1-q)/4),	((1-q)/4),	0.322*q, 	0.298*q, 	0.078*q, 	0.302*q   ;
    	((1-q)/4),	((1-q)/4),	((1-q)/4),	((1-q)/4),	0.248*q, 	0.246*q, 	0.298*q, 	0.208*q   ;
    	((1-q)/4),	((1-q)/4),	((1-q)/4),	((1-q)/4),	0.177*q, 	0.239*q, 	0.292*q, 	0.292*q   ;  ];

%emission matrix
%       |A+|C+|G+|T+|A-|C-|G-|T-|
EMIS = [ 1, 0, 0, 0, 1, 0, 0, 0 ; 
         0, 1, 0, 0, 0, 1, 0, 0 ;
         0, 0, 1, 0, 0, 0, 1, 0 ;
         0, 0, 0, 1, 0, 0, 0, 1 ; ];
    
%default variables    
cpgIslands = [];
cpgCount = 0;
cpgLength = 0;
cpgLengthMtr = [];
mtrCount = 1;
codingCount = 0;

%%%%%%%%5%%%%%%%%%%%%%%
%the viterbi algorithm
%%%%%%%%%%%%%%%%%%%%%%%

states = 8;
%seqx = toy1;
len = length(seqx);
% take the logs of the matrices
logE = log(EMIS);
logT = log(TRANS);

% allocate space
ptr_v = zeros(states,len);
% initialize the matrix
% the model is in state 1 at step 0
start_prob=log((1/states));
v = zeros(states,1);
v=v+start_prob;
%v(1,1) = 0;
placeholder = v;

% iteration step
for i = 1:len
    for st = 1:states
        % Vj(i) = ej(xi) × maxk akj Vk(i-1) 
        maxValue = -inf;
        maxPtr = 0;
        % keeps track of max while looping
        for in = 1:states 
            val = placeholder(in) + logT(in,st);
            if val > maxValue
                maxValue = val;
                maxPtr = in;
            end
        end
        %termination
        %P(x, π*) = maxk Vk(N) 
        ptr_v(st,i) = maxPtr;
        
        %indexing g,c,t,a
        if (seqx(i) == 'G')
            pos = 3;
        elseif (seqx(i) == 'T')
            pos = 4;
        elseif (seqx(i) == 'C')
            pos = 2;
        elseif (seqx(i) == 'A')
            pos = 1;
        % ignore if unknown nucleotide
        elseif (seqx(i) == 'N')
           % pos = 1;
        end
        % update the matrix
        v(st) = logE(pos,st) + maxValue;
    end
    placeholder = v;
end

% traceback step initialization
[tempVal, final] = max(v);
current = zeros(1,len);
bp = 0;
trackEnd = length(seqx);
trackStart = 0;
current(len) = final;

%~3.4M basepairs from 43507093 to 4694432
chr21Start = 43507093;
%scan the known coding regions
fid = fopen('Chr21.txt', 'r');
formatSpec = '%d%d%c%s';
genes = textscan(fid,formatSpec, 'delimiter', '\n');
gene1length = length(genes{1});

%scan the known motifs
%fid2 = fopen('motifs_Chr21.txt', 'r');
fid2 = fopen('temp_motifs.txt', 'r');
formatSpec2 = '%s%s%d%d%c';
motifs = textscan(fid2,formatSpec2, 'delimiter', ' ');
motif1length = length(motifs{1});

%start of traceback step
for count = len-1:-1:1
    current(count) = ptr_v(current(count+1),count+1);
    
    %identifying cpg islands (in the + region)
    if (current(count) < 5 && (count > 2))
        bp = bp + 1;
    %if more than 200 nucleotides in a row, declare it a CpG island
    elseif (bp >= 200)
        cpgCount = cpgCount + 1;
        trackStart = count;
        cpgLength = trackEnd - trackStart;
        %keep track of the CpG island start/end points for motif detection
        cpgLengthMtr(mtrCount, 1) = trackStart;
        cpgLengthMtr(mtrCount, 2) = trackEnd;
        cpgLengthMtr(mtrCount, 3) = cpgCount;
        mtrCount = mtrCount + 1;
        
        %display as text string -- displays from last island to first because of backwards traceback
        cpgIslands = [' CpG Island ', num2str(cpgCount), ': ', num2str(cpgLength), 'bp (', num2str(trackStart), ' - ', num2str(trackEnd), ') '];

        %for loop here for gene identification (within 500bps)
        for k = 1:500
            %initialize the positions to account for chr21 starting
            %position
            pos_position = (trackStart - k) + chr21Start;
            neg_position = (trackEnd + k) + chr21Start;
            
            %genes{1}(1) == start of gene | genes{4}{1} == gene names
            %for all the first indices in genes
            for posit = 1:gene1length
                %if the position matches the start of a gene (identifying in reverse direction)
                if (pos_position == genes{1}(posit))
                    codingCount = codingCount + 1;
                    %concatentate the output to include the gene
                    cpgIslands = [cpgIslands genes{4}{posit} ' (' num2str(genes{3}(posit)) ')'];
                    
                
                %identifying genes in the forward direction
                elseif (neg_position == genes{1}(posit))
                    codingCount = codingCount + 1;
                    cpgIslands = [cpgIslands genes{4}{posit} ' (' num2str(genes{3}(posit)) ')'];
                end
            end
            %{
            % motif detection loop
            bp_start = (trackStart + k) + chr21Start;
            %bp_end   = cpgLengthMtr(k, 2) + chr21Start;
            tempCnt = 0;
            for posit2 = 1:motif1length              
                if ((bp_start == motifs{3}(posit2))) %&& (motifs{4}(posit) <= bp_end))
                    %the CpG island it's mapped to
                    motif_out = [motifs{1}{posit2} '; CpG Island ' num2str(cpgCount)];
                    %displays each individual motif along with the island
                    disp(motif_out);
                    tempCnt = tempCnt + 1;
                end
            end
            %displays the hits on each CpG island
            %cpgIslands = [cpgIslands '; ' num2str(tempCnt)];
            %}
        end
        %output the Islands one-by-one
        disp(cpgIslands); 
        bp = 0;
        trackEnd = 0;
    else
        bp = 0;
        %starts at end
        trackEnd = count;
    end
end
       
% Proper output format
% Total CpG islands found: 3 ; 2 out of 3 islands are followed by a coding region
out1 = 'Total CpG islands found: ';
out2 = ' out of ';
out3 = ' islands are followed by a coding region';
output1 = [out1, num2str(cpgCount), ' ; ', num2str(codingCount), out2, num2str(cpgCount), out3];
disp(output1);
fclose(fid);
fclose(fid2);

