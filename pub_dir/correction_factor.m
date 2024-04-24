function out = correction_factor(win,N2)
% function out = correction_factor(win,N2)
% Gives the correction factor for the estimator described in [1,2,3] 
% for an arbitrary window 
%
% [1] C. Candan, "A Method For Fine Resolution Frequency Estimation From Three
% DFT Samples," IEEE Signal Processing Letters, Vol. 18, No.6, p. 351-354,
% June 2011. 
% 
% [2] C. Candan, "Analysis and Further Improvement of Fine Resolution 
% Frequency Estimation Method From Three DFT Samples," 
% IEEE Signal Processing Letters, vol.20, no.9, pp.913–916, Sept. 2013.
%
% [3] C. Candan, "Fine Resolution Frequency Estimation From Three DFT
% Samples: Windowed Case," Elsevier Signal Processing, vol. 114, p. 245-250, 
% Sept. 2015.
%
%
% Input: 
% win : Windowing function 
% N2  : Number of DFT points 
%
% Output: 
% out : Correction factor
%
% Note: 
% -----
% The case of zero-padding (N2 > length(win)) is identical to the
% correction factor of the zero-padded window, i.e. 
% win_zero_padded = [win(:) zeros(N2-length(win),1)]; 
% 
% Some usage examples: 
% win = hamming(N)';
% win = window(@hamming,N)';
% win = window(@blackman,N)';
%
% c_N = correction_factor(hamming(16),32); 
%
% Oct. 2013, 
% Cagatay Candan
%

win = win(:).';

if exist('N2')==0, 
    N2 = length(win); 
elseif N2~=length(win),
    win = [win zeros(1,N2-length(win))];
end;

N = length(win); 
nvec = 0:N-1;

fd = @(inp) sum(repmat(win,length(inp),1).*exp(j*2*pi/N*inp(:)*nvec),2); 
fdp = @(inp) j*2*pi/N*sum(repmat(win.*nvec,length(inp),1).*exp(j*2*pi/N*inp(:)*nvec),2); 

A0 = imag(fd(1)-fd(-1));
A1 = fdp(1) - fdp(-1);

B0 = 2*fd(0) - fd(1) - fd(-1); 
B1 = imag(2*fdp(0) - fdp(1) - fdp(-1)); 
%%%%%%%

c1 = (A1*B0+A0*B1)/B0^2;
out = 1/c1;

