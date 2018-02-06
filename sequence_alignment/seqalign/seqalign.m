%humanHEXA = getgenbank('NM_000520.5');
%mouseHEXA = getgenbank('AK080777');

humanHEXA= fastaread('NM_000520.fasta'); 
mouseHEXA= fastaread('AK080777.fasta'); 

%load hexosaminidase;
humanORFs = seqshoworfs(humanHEXA.Sequence);
mouseORFs = seqshoworfs(mouseHEXA.Sequence);

humanProtein = nt2aa(humanHEXA.Sequence);
mouseProtein = nt2aa(mouseHEXA.Sequence);

%dot plot
seqdotplot(mouseProtein,humanProtein);
xlabel('Human hexosaminidase A');ylabel('Mouse hexosaminidase A');

seqdotplot(mouseProtein,humanProtein,10,8)
xlabel('Human hexosaminidase A');ylabel('Mouse hexosaminidase A');

[score, globalAlignment] = nwalign(humanProtein,mouseProtein);
showalignment(globalAlignment);

[score, globalAlignment] = nwalign(humanProtein,mouseProtein,'glocal',true);
showalignment(globalAlignment);

%translate to proteins
humanStart = find(humanProtein == 'M',1);
mouseStart = find(mouseProtein == 'M',1);
humanStop = find(humanProtein(humanStart:end)=='*',1) + humanStart - 1;
mouseStop = find(mouseProtein(mouseStart:end)=='*',1) + mouseStart - 1;
humanSeq = humanProtein(humanStart:humanStop);
humanSeqFormatted = seqdisp(humanSeq)
mouseSeq = mouseProtein(mouseStart:mouseStop);
mouseSeqFormatted = seqdisp(mouseSeq)
%%%%%%%%%%%

%global alignment of proteins
[score, alignment] = nwalign(humanSeq,mouseSeq);
showalignment(alignment);

%local alignment of proteins
[score, localAlignment] = swalign(humanProtein,mouseProtein);
showalignment(localAlignment);





