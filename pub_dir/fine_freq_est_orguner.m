function outfreq = fine_freq_est_orguner(data,maxorder);
% function outfreq = fine_freq_est_orguner(data,maxorder);
% Generates frequency estimates via fusing maxorder estimates given by 
%
% \deltahat_k = N/pi atan(tan(pi k/N) Real { Ratio } )
%
% where 
% Ratio = (R[kp+k]e^{-jk pi/N} - R[kp - k]e^{+jk pi/N}) /
%        (R[kp+k]e^{-jk pi/N} + R[kp - k]e^{+jk pi/N} - 2R[kp]/cos(pi k/N))
%
% data :  N x Mcnum matrix ,(N: number of samples, Mcnum: number of vectors)
% maxorder : Number of estimates to be fused (maxorder < N/2 )
% outfreq : 1 x Mcnum vector
%
% March 2014
% 

[N,MCnum] = size(data);
if exist('maxorder')==0, maxorder = N/2 - 1; end; 

umutall = zeros(maxorder,MCnum);
%%%First Estimate
morder = 1; %do the first estimation (\deltahat_1)
[outfreq_dum,outf,outmaxind] = est_sub(data,morder);
umutall(morder,:) = outfreq_dum;
clear data; 

%%%Remaining Estimates
for morder=2:maxorder,  %Do the rest, if maxorder>1
    outfreq_dum = est_sub([],morder,outf,outmaxind);
    umutall(morder,:) = outfreq_dum;
end;

%Fusion
dum = 1:maxorder; fusw = 1./(sin(pi*dum/N)).^2; fusw = fusw/sum(fusw); 
outfreq = fusw*umutall;


%%%%%%%%%%%%%%%
function [outfreq,outf,outmaxind] = est_sub(data,morder,outf,outmaxind);

if exist('outf')==0, %If outf is provided at the input, do not repeat 
                     %                                FFT calculation
    [N,MCnum] = size(data);
    outf=fft(data,[],1); clear data;
    out = real(outf).*real(outf) + imag(outf).*imag(outf);
    [outmaxval,outmaxind] = max(out,[],1);
else
    [N,MCnum] = size(outf);
end;

dumvec = (0:length(outmaxind)-1)*N; 
vec = outmaxind -   1;  % 0 <= vec elements <= N-1
vecp = mod(vec + morder, N); % 0 <= vec elements <= N-1 
vecm = mod(vec - morder, N); % 0 <= vec elements <= N-1 

vec = vec + dumvec + 1; 
vecp = vecp + dumvec + 1;
vecm = vecm + dumvec + 1;
index = [vecm;vec;vecp];
out  = outf(index); 

%
%Constants
expp=exp(j*pi*morder/N); expm=exp(-j*pi*morder/N);
mycos=cos(pi*morder/N); mytan=tan(pi*morder/N); 

%
out = [-expp 0 expm; expp -2/mycos expm]*out; 
outfreq = N/pi*atan(real(mytan*out(1,:)./out(2,:)));
outfreq = outfreq + outmaxind - 1;