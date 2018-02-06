FASTAfilename = 'seq1.out';
GenesFilename = 'Chr21.txt';
MotifsFilename = 'motifs_Chr21.txt';

min_size = 200;

if ~isempty(strfind(FASTAfilename, '.out'))
    fidIn = fopen(FASTAfilename, 'r');
    seq = textscan(fidIn, '%s');
    seq = strjoin(seq{1}','');
    fclose(fidIn);
elseif ~isempty(strfind(FASTAfilename, '.fasta'))
    seq = fastaread(FASTAfilename);
    seq = seq.Sequence;
end

% define p and q (you might need to tune these values)
p = 0.98;
q = 0.999;
disp('Performing Viterbi...');
CpGs = viterbi(seq, p, q); 
disp('Finished Viterbi');

% read genes
if ~isempty(GenesFilename)
    genes = readtable(GenesFilename);
    genes.Properties.VariableNames = {'Start' 'End' 'Direction' 'Name'};
end

% read motifs
if ~isempty(MotifsFilename)
    motifs = readtable(MotifsFilename);
    motifs.Properties.VariableNames = {'Name' 'Chr' 'Start' 'End' 'Direction'};
end

for i=1:size(CpGs,1)
    
    if (CpGs(i,2)-CpGs(i,1)+1) >= min_size
        % CpG size greater than min_size: 
        % - report CpG
        
        
        % - Identify the down-stream genes that immediately follow your 
        % CpG islands (starting within 500 basepairs).
        
        
        % - Map the motifs to your CpG islands and report the motif hits.
        
    end
    
end




