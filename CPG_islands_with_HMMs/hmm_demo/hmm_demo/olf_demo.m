orseqs = fastaread('347OR.fasta'); 

or1 = orseqs(1).Sequence

%proteinplot(or1)

open('hydrophobicity.fig'); % Comment this if you want to use the proteinplot

hold on

for i=1:20
    intseqs(i) = {aa2int(orseqs(i).Sequence)};
end

%T = [0.6 0.4;
%     0.4 0.6];


T = [0.8 0.2;
     0.2 0.8];

E = [0.018	0.067	0.067	0.067	0.018	0.067	0.067	0.067	0.067	0.01	0.01	0.067	0.018	0.018	0.067	0.067	0.067	0.067	0.067	0.01;
     0.114	0.007	0.007	0.007	0.114	0.007	0.007	0.025	0.007	0.114	0.114	0.007	0.114	0.114	0.025	0.025	0.025	0.025	0.025	0.114];

[estT, estE] = hmmtrain(intseqs, T, E, 'maxiterations',40)

estimatedStates = hmmviterbi(aa2int(or1),estT,estE);
plot(estimatedStates)
hold off