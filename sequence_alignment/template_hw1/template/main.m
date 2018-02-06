%% Delete sequence_analysis.txt

%this file will be created later on to report sequence alignment scores,
%and alignments but it will be deleted in order to prevent un-needed
%appending

%delete('sequence_analysis.txt')

%% load sequences and matched sequences

seq1 = fastaread('Human_HOX.fa');
seq2 = fastaread('Fly_HOX.fa');

%to label which alignment the score and alignment come from 
proteinAligned='HOX Alignment';
%matched regions for anchored needlman wunsch later
matched_regions = dlmread('Match_HOX.txt');

seq1 = seq1.Sequence;

seq2 = seq2.Sequence;


%% penalties
gap = -2;
match = 1;
mismatch = -3;

%% needleman wunsch

 [score_need_wunsch, alignment1_need_wunsch, aligner_need_wunsch, alignment2_need_wunsch]...
     = needleman_wunsch(seq1, seq2, match, mismatch, gap,matched_regions);

 transpose_need_wunsch_Alignment=[alignment1_need_wunsch';aligner_need_wunsch';alignment2_need_wunsch'];
 
% alignmentType='needleman wunsch';
 
%  dlmwrite('sequence_analysis.txt',proteinAligned,'newline','unix','-append','delimiter','')
%  dlmwrite('sequence_analysis.txt',alignmentType,'-append','newline','unix','delimiter','') 
%  dlmwrite('sequence_analysis.txt',score_need_wunsch,'newline','unix','-append','delimiter','')
%  dlmwrite('sequence_analysis.txt',transpose_need_wunsch_Alignment,'-append','newline','unix','delimiter','') 
%  

%% anchored needleman wunsch
[score_anc_need_wunsch, alignment1_anc_need_wunsch, aligner_anc_need_wunsch, alignment2_anc_need_wunsch]...
    = anchored_needleman_wunsch(seq1, seq2, match, mismatch, gap,matched_regions);

transpose_anc_need_wunsch_Alignment=[alignment1_anc_need_wunsch';aligner_anc_need_wunsch';alignment2_anc_need_wunsch'];
 
%alignmentType='anchored needleman wunsch';
 
% dlmwrite('sequence_analysis.txt',proteinAligned,'newline','unix','-append','delimiter','')
%  dlmwrite('sequence_analysis.txt',alignmentType,'-append','newline','unix','delimiter','') 
%  dlmwrite('sequence_analysis.txt',score_anc_need_wunsch,'newline','unix','-append','delimiter','')
%  dlmwrite('sequence_analysis.txt',transpose_anc_need_wunsch_Alignment,'-append','newline','unix','delimiter','') 

%% perform alignment of random sequences 10,000 times
% acc_scores = [];
% for i=1:1e4  
% %for i=1:100
%  % permute seq1 and create loop in order to determine new amino
%     % acid positions
%      permute_seq1 = randperm(length(seq1));
%      %create character array for new sequence
%      seq1_rand= [''];
%      for aa_pos=1:length(seq1)
%          seq1_rand(aa_pos)=seq1(permute_seq1(aa_pos));
%      end
%      
%      %and repeat for seq2 
%      permute_seq2=randperm(length(seq2));
%      seq2_rand = [''];
%      for aa_pos=1:length(seq2)
%          seq2_rand(aa_pos)=seq2(permute_seq2(aa_pos));
%      end
%      
%      %now once finished seq1 = seq1_rand and seq2 = seq2_rand
%     seq1=seq1_rand;
%     seq2=seq2_rand;
%      
%     % call anchored_needleman_wunsch and store scores into acc_scores
%     [permute_score_anc_need_wunsch, permute_alignment1_anc_need_wunsch, permute_aligner_anc_need_wunsch, permute_alignment2_anc_need_wunsch]= anchored_needleman_wunsch(seq1, seq2, match, mismatch, gap,matched_regions);
%     acc_scores=[acc_scores;permute_score_anc_need_wunsch];
% end
% % 
% %% plot random scores and actual score (uncomment below to plot)
% figure1 = figure;
% axes1 = axes('Parent',figure1);
% histogram(acc_scores)
% hold on
% %histogram(score)
% histogram(score_anc_need_wunsch)
% grid on; grid minor;
% legend({'Random Sequences', 'Alignment Score'});
% xlabel('Score'); ylabel('Frequency');
% set(axes1,'FontSize',14);
