Submitted by Gautam Shenoy. Dont Spam me. 


Description:
This program describes the use of Viterbi algorithm for deducing the most probable path given a set of observations
and to deduce hidden states. Hope you find it useful.



Notes:
1) dec_state will consist of a set of numbers 1=> Sunny, 2=> Cloudy and 3=> Rainy

2) decode_state gives the final o/p in words

3) Got the algorithm from Rabiner's paper entitled "A tutorial on HMM and
 application in speech recognition".

4) This line is not necessary. But I wanted to see if the results
 would be different If I attempted a decoding without backtracking i.e. as I progress. For
 what I tried the results were not always same. This example shows that mx is not the same as
dec_state.


5) The TPM(Transition Probability matrix) [T_prob] should have a rowsum of 1. i.e. 
for every row the sum of the columns should equal 1. But this one has a columnsum of 1. 
But that is the matrix given on the site. It should have had a rowsum of 1. 
But If I use the transpose of this the answer comes out wrong. So for verification purpose 
I used their matrix and the answer came out right.

6) Example courtesy of www.comp.leeds.ac.uk. (Just type viterbi algorithm on Google and youll find it).

For your application, use a rowsum=1 matrix for tpm(AKA stochastic matrix).
        
Keywords: Viterbi, Hidden markov model      
        
Any discrepancies contact me at gau_dmc_engg@yahoo.co.in