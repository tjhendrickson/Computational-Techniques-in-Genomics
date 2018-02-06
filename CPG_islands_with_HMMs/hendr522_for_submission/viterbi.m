function CpGs =viterbi(seq, p, q)  
% viterbi Perform the Viterbi algorithm.
%
% Input variables:
% seq: string containing the sequence
% p: probability to continue in the CpG island
% q: probability to continue out of the CpG island
%
% Output variables:
% CpGs: contains a list of pairs identifying CPG islands
%
% Output example:
% [20    412    <--- CpG1 (from position 20 to 412 including)
%  1027  1456   <--- CpG2 (from position 1027 to 1456 including)
%  2213  2560]  <--- CpG3 (from position 2213 to 2560 including)
%
%%%%%%%%%%%%%% YOUR CODE STARTS HERE

format short

%insert transition probabilities matrix from homework2 webpage  
trans_prob=[0.180*p, 0.274*p, 0.426*p, 0.120*p, ((1-p)/4), ((1-p)/4),((1-p)/4),((1-p)/4);
            0.171*p, 0.368*p, 0.274*p, 0.188*p, ((1-p)/4), ((1-p)/4), ((1-p)/4), ((1-p)/4);
            0.161*p, 0.339*p, 0.375*p, 0.125*p, ((1-p)/4), ((1-p)/4), ((1-p)/4), ((1-p)/4);
            0.079*p, 0.355*p, 0.384*p, 0.182*p, ((1-p)/4), ((1-p)/4), ((1-p)/4), ((1-p)/4);
            (1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4,	0.300*q, 0.205*q, 0.285*q, 0.210*q;
            (1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4, 0.322*q, 0.298*q, 0.078*q, 0.302*q;
            (1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4,	0.248*q, 0.246*q, 0.298*q, 0.208*q;
            (1-q)/4, (1-q)/4, (1-q)/4, (1-q)/4,	0.177*q, 0.239*q, 0.292*q, 0.292*q];

%now insert emission probabilities for all states. The order is as follows:
%first four rows correspond to A,C,G, T within CG rich area, and second four
%rows correspond to A,C,G,T in non CG rich area

emis_prob=[ 1, 0, 0, 0, 1, 0, 0, 0 ; 
         0, 1, 0, 0, 0, 1, 0, 0 ;
         0, 0, 1, 0, 0, 0, 1, 0 ;
         0, 0, 0, 1, 0, 0, 0, 1 ; ]; 
%how many states are there?
states=8;

%initialize variables
start_prob=log((1/states));
score=zeros(states,1);
score=score+start_prob;
temp=score;


% create a variable of zeroes which corresponds to the number of states as
% rows and number of observations as columns 
V=zeros(states,(length(seq)));

%now create reference row corresponding to sequence location

V(1,1:end)=(1:length(seq));

for i=1:length(V)
            %determine the nucleotide and it's corresponding emission position
        if strcmp(seq(V(1,i)),'A')
           pos=1;

        elseif  strcmp(seq(V(1,i)),'C')
            pos=2;


        elseif strcmp(seq(V(1,i)),'G')
            pos=3;


        elseif strcmp(seq(V(1,i)),'T')
            pos=4;

        end
        %set emission row equal to pos 1:4 corresponding to its correct
        %nucleotide
        emission=emis_prob(pos,:);
        holder=zeros(states,1);
        for j=1:states
                max_ptr=0;
                for st=1:states
                    value(st,j)=((log(trans_prob(st,j)))) + temp(st);
                    %if value > maxVal
                    %    maxVal=value;
                    %    max_ptr=st;
                    %end
                
                end
                % determine score for state j
                score(j)=max(value(:,j))+log(emission(j));
                %determine value and pointer for particular state j
                [maxVal,max_ptr]=max(value(:,j));
                
                %place pointer into variable holder
                holder(j,:)=max_ptr;
        end
        %set holder equal to V(i) matrix
        V(:,i)=holder;
        %set temp to score
        temp=score;
end


%% now perform traceback to find CpG islands

%what is the max position and value of final column? 
[value,final_position]=max(score);

%initialize variables
basePairs=0;
Ptr=zeros(1,length(V));
Ptr(length(V))=final_position;
CpG_end=length(V);
CpG_start=0;
%set dummy variable to count CpG islands
x=1;
CpGs=[];
for traceback = length(V)-1:-1:1
    Ptr(traceback)=V(Ptr(traceback+1),traceback+1);
   if Ptr(traceback) < 5 % we are in a CpG island
       %add to base pairs
       basePairs=basePairs+1;
   else
       
       if basePairs > 3
           %now that CpG island has ended determine start point
           CpG_start=traceback;

           CpGs(x,:)=[CpG_start,CpG_end];
           
           x=x+1;
       end

       %no longer in CpG island so set base pair count to zero
       basePairs=0;
       CpG_end=traceback;
       
       
   end
end
end