clear

%FASTAfilename = 'seq1.out';
FASTAfilename = 'HMC21_NT_011515.fasta';
GenesFilename = 'Chr21.txt';
MotifsFilename = 'motifs_Chr21.txt';

outputFile='output_cpgIslands.txt';
if exist(outputFile)
    delete(outputFile)
end

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
p = 0.99;
q = 0.98;
disp('Performing Viterbi...');
CpGs =[];
CpGs =viterbi(seq, p, q);   
disp('Finished Viterbi');
%% Now create CpG output

% read genes
if ~isempty(GenesFilename)
    genes = readtable(GenesFilename);
    genes.Properties.VariableNames = {'Start' 'End' 'Direction' 'Name'};
end

% read motifs
if ~isempty(MotifsFilename)
    motifs = readtable(MotifsFilename);
    motifs.Properties.VariableNames = {'Name' 'Chr' 'Start' 'End' 'Direction'};
    motifs.Direction=categorical(motifs.Direction);
end
count=1;
coding_regions=0;
chr21_start=43507093;
chr21_end=46944323;
islandHolder=cell([(length(CpGs)),1]);
for i=1:size(CpGs,1)
    
    
    if (CpGs(i,2)-CpGs(i,1)) >= min_size
        
        
        % CpG size greater than min_size: 
        
        %how many base pairs long is the CpG island?
        basePairs=(CpGs(i,2)-CpGs(i,1));
        
        %start and stop CpG island locations
        start_island=CpGs(i,1)+chr21_start;
        end_island=CpGs(i,2)+chr21_start;
        
        %creating output to report CpG island
        
        island=['CpG island ', num2str(count), ': ', num2str(basePairs), ' bp (', num2str(start_island) '-', num2str(end_island), ')'];
        
        %% - Identify the down-stream genes that immediately follow your CpG islands (starting within 500 basepairs).
        
        % ending point of gene hunt downstream from CpG island
        geneHunt=end_island+500;
        negStrand_chr21=geneHunt;
        posStrand_chr21=end_island;
        for bp=1:500
            for gene=1:height(genes)
                 if strcmp(genes.Direction(gene),'+')                 
                        genehit=posStrand_chr21==genes.Start(gene);
                        if genehit==1
                            name=char(genes.Name(gene));
                            comment='; gene name';
                            geneReport=[comment ' ' name];

                            %append to island variable
                            island=[island ' ' geneReport];

                            coding_regions=coding_regions+1;
                        end
                 elseif  strcmp(genes.Direction(gene),'-')
                    genehit=negStrand_chr21==genes.End(gene);
                    if genehit==1
                        name=char(genes.Name(gene));
                        comment='; gene name';
                        geneReport=[comment ' ' name];

                        %append to island variable
                        island=[island geneReport];

                        coding_regions=coding_regions+1;
                    end                    
                 end
            end 
            posStrand_chr21=(bp+end_island);
            negStrand_chr21=(geneHunt-bp);
        end
        
        %% motif analysis
        %if motif start and end is within cpg island return a 1 for that
        %motif
        motif=motifs.Start>=CpGs(i,1)+chr21_start & motifs.End<=CpGs(i,2)+chr21_start;
        vars={'Name'};
        %create column with motif hit names
        motifhits=motifs(motif,vars);
        
        %loop through motif hits and append to CpG island text
         for x=1:height(motifhits)

                    motifName=char(motifhits.Name(x));
                    comment=' ; motif name';
                    motifReport=[comment ' ' motifName];
                    island=[island ' ' motifReport];
         end

        
        %convert from string to cell
        island=cellstr(island);
        %put CpG island into cell array islandHolder
        islandHolder(count)=island;
        count=count+1;
    end       
end
%% Output
%reduce cell array so that only data is in each cell
 island_all=islandHolder(~cellfun('isempty',islandHolder));
 
%Total CpG islands found: 3 ; 2 out of 3 islands are followed by a coding region  
report_cpgs=[ 'Total CpG islands found: ', num2str(count-1), ' ; ', num2str(coding_regions), ...
    ' out of ', num2str(count-1), ' islands are followed by a coding region.'];

%final output
final_output=[report_cpgs; island_all];

%convert final output to character array for easier reading
final_output=char(final_output);

%now write out put to a text file
dlmwrite(outputFile,final_output,'delimiter', '')





