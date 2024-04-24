function [outfreq2,outfreq1] = windowed_fine_freq_est(data,win,Npoint);
%function [outfreq2,outfreq1] = windowed_fine_freq_est(data,win,Npoint);
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
% [3] C. Candan, "Fine Resolution Frequency Estimation From Three DFT
% Samples: Windowed Case," Elsevier Signal Processing, vol. 114, p. 245-250, 
% Sept. 2015.
%
% data :  N x Mcnum matrix ,(N: number of samples, Mcnum: number of vectors)
% win  :  window to be used 
%         win = hamming(N);
%         win = window(@hamming,N);
%         win = window(@blackman,N);
%         win = window(@rectwin,N);
%         win = window(@hann,N);
%         win = window(@gausswin,N,2.5); 
%
% Npoint : Number of DFT points (shown as N2 in [3]) 
% outfreq2 : 1 x Mcnum vector (outfreq is in FFT bins) 
%            uses \hat{\delta} = \hat{\delta}_1 + \hat{\delta}_2, 
%            (See Table 1 and Figure 3 of [3]) 
% outfreq1 : 1 x Mcnum vector (outfreq is in FFT bins)
%            uses \hat{\delta} = \hat{\delta}_1, 
%
%
% C. Candan
% Oct. 2013
%

[N,MCnum] = size(data);
if length(win)~=N, 
    disp('Window length should match input length (N)'); 
return; end;

if exist('Npoint')==0, 
    Npoint = N; 
elseif Npoint~=N,
    data=[data; zeros(Npoint-N,MCnum)];
    win = [win(:); zeros(Npoint-N,1)];
end;

cf = correction_factor(win);
[deltahat1,outmaxind] = est_delta(data,win,cf,Npoint,MCnum); %INITIAL ESTIMATE

data = data.*exp(-j*2*pi/Npoint*(0:(size(data,1)-1))'*deltahat1);
[deltahat2,outmaxind] = est_delta(data,win,cf,Npoint,MCnum); %SECOND ESTIMATE

outfreq1  = deltahat1+ outmaxind - 1;  
outfreq1  = outfreq1*N/Npoint; 

outfreq2  = deltahat1+ deltahat2 + outmaxind - 1;  %COMBINED ESTIMATE
outfreq2  = outfreq2*N/Npoint; 



%%%%%
function  [deltahat,outmaxind] = est_delta(data,win,cf,N,MCnum);

outf=fft(repmat(win(:),1,MCnum).*data,[],1); 
out = real(outf).*real(outf) + imag(outf).*imag(outf);
[outmaxval,outmaxind] = max(out,[],1); 
dumvec = (0:length(outmaxind)-1)*N; 
vec   = outmaxind + dumvec;
vecp1 = outmaxind + 1; vecp1(vecp1==(N+1))=1; vecp1 = vecp1 + dumvec;
vecm1 = outmaxind - 1; vecm1(vecm1==0)=N; vecm1 = vecm1 + dumvec;
index = [vecm1;vec;vecp1]; clear dumvec

out  = outf(index); clear outf; 
%
out  = [1 0 -1; -1 2 -1]*out; 
outfreq = real(out(1,:)./out(2,:)); 
deltahat = cf*outfreq; 

%%%end of function est_delta


%%%%%%%%%%%%%%%%%%%
function cf = correction_factor(win)
% Gives the correction factor for Jacobsen type estimator for an arbitrary
% window 
%
% win : window function 
% 
% some examples:
% win = hamming(N)';
% win = window(@hamming,N)';
% win = window(@blackman,N)';
% win = window(@rectwin,N)';
% win = window(@hann,N)';
% win = window(@gausswin,N,2.5); 
%
%Sept. 2013, 
%Cagatay Candan
%

win = win(:)'; N = length(win); nvec = 0:N-1;

fd = @(inp) sum(repmat(win,length(inp),1).*exp(j*2*pi/N*inp(:)*nvec),2); 
fdp = @(inp) j*2*pi/N*sum(repmat(win.*nvec,length(inp),1).*exp(j*2*pi/N*inp(:)*nvec),2); 
%fdpp = @(inp) (j*2*pi/N)^2*sum(repmat(win.*nvec.*nvec,length(inp),1).*exp(j*2*pi/N*inp(:)*nvec),2); 

A0 = imag(fd(1)-fd(-1));
A1 = fdp(1) - fdp(-1);
%A2 = imag((fdpp(1) - fdpp(-1))/2);

B0 = 2*fd(0) - fd(1) - fd(-1); 
B1 = imag(2*fdp(0) - fdp(1) - fdp(-1)); 
%B2 = (2*fdpp(0) - fdpp(1) - fdpp(-1))/2;

%%%%%%%
c1 = (A1*B0+A0*B1)/B0^2;
%c3 = (A1*B2+A2*B1)/B0^2;

cf = 1/c1;