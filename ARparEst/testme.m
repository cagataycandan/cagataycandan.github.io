A = [1 -2.7607 3.8106 -2.6535 0.9238];  %
  % Filter coefficients are from
  % S. Kay, Recursive maximum likelihood estimation of autoregressive processes,
  %  IEEE Transactions on Acoustics, Speech, and Signal Processing
  %  Vol.31, p. 56-65, 1983. 
b = 0.36; 
N = 50; %Observation vector length

ARorder = length(A) - 1; %No model mismatch

y = filter(1,A,b*randn(N,1)); %generate observation vector
[a_wFB,b_wFB,a_FB,b_FB] = wFBprediction(y,ARorder); %1st stage output
[anew,bnew] = approx_MLc(y,a_wFB); %2nd stage output with a_WFB as initial condition

%display loglikelihoods
str1 = sprintf('Loglikelihoods: \n \t FB prediction (conventional)\t :%0.7g',loglikelihood(a_FB,b_FB,y));
str2 = sprintf('\n \t 1st stage output (weighted FB)\t :%0.7g',loglikelihood(a_wFB,b_wFB,y));
str3 = sprintf('\n \t 2nd stage output\t\t\t\t :%0.7g',loglikelihood(anew,bnew,y));
sprintf('%s %s %s',str1, str2, str3),