function [outfreq_candan_rem, outfreq_candan_cor,outfreq_jac] = fine_freq_est(data);
%function outfreq = fine_freq_est(data);
% Implements fine frequency estimation method described in 
%
% [1] C. Candan, "A Method For Fine Resolution Frequency Estimation From Three
% DFT Samples," IEEE Signal Processing Letters, Vol. 18, No.6, p. 351-354,
% June 2011. 
% 
% [2] C. Candan, "Analysis and Further Improvement of Fine Resolution 
% "Frequency Estimation Method From Three DFT Samples," 
% IEEE Signal Processing Letters, vol.20, no.9, pp.913–916, Sept. 2013. 
%
% data :  N x Mcnum matrix ,(N: number of samples, Mcnum: number of vectors)
% outfreq : 1 x Mcnum vector (outfreq is in FFT bins)
%
% outfreq_candan_rem : bias-removed estimate ([2])
% outfreq_candan_cor : bias-corrected estimate ([1,2])
% outfreq_jac : Jacobsen's estimator ([1])
%
% C. Candan
% Sept. 2013
%

[N,MCnum] = size(data);

outf=fft(data,[],1); clear data; 
out = real(outf).*real(outf) + imag(outf).*imag(outf);
[outmaxval,outmaxind] = max(out,[],1); 
dumvec = (0:length(outmaxind)-1)*N; 
vec   = outmaxind + dumvec;
vecp1 = outmaxind + 1; vecp1(vecp1==(N+1))=1; vecp1 = vecp1 + dumvec;
vecm1 = outmaxind - 1; vecm1(vecm1==0)=N; vecm1 = vecm1 + dumvec;
index = [vecm1;vec;vecp1]; clear dumvec

out  = outf(index); clear outf

%
out  = [1 0 -1; -1 2 -1]*out; 
outfreq_jac = real(out(1,:)./out(2,:)); 
outfreq_candan_cor = tan(pi/N)/(pi/N)*outfreq_jac; 
outfreq_candan_rem = atan(outfreq_candan_cor*pi/N)*N/pi;

%
outfreq_jac = outfreq_jac + outmaxind - 1; 
outfreq_candan_cor = outfreq_candan_cor + outmaxind - 1; 
outfreq_candan_rem = outfreq_candan_rem  + outmaxind - 1; 
