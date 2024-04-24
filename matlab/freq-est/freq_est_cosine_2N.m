function deltahat_fused = freq_est_cosine_2N(rmat,maxiter)
%function est_fused = freq_est_cosine_2N(rmat,maxiter)
%
% Returns the frequency estimate for a real-valued sinusoidal observed under
% AWGN. 
% 
% Implements the method given in 
%      "A Computationally Efficient Fine Frequency Estimation Method 
%                      For Real-Valued Sinusoids"
% by Cagatay Candan and Utku Celebi.
%   
% Inputs:
%     rmat:    N x MCrun matrix where each column contains an observation
%              vector of N samples
%     maxiter: maximum number of iterations for the inverse mapping of the
%              invariant function, (default value, if not provided, = 6)
%Output: 
%     est_fused : The frequency with the unit of N-point DFT bins, i.e. 
%                 a real number in [0,N/2] where N is number of samples 
%                 for each observation.
%                   
%                 To convert est_fused to radian per samples, use
%                             omega  = 2*pi/N*est_fused; 
%
%Sept. 2018
%CC
%
if nargin==1, maxiter = 6; end; 

[N,trial] = size(rmat);

R =fft(rmat,2*N,1);
[~,max_frequency_bins]=max(real(R(1:N,:)).^2 + imag(R(1:N,:)).^2,[],1);  %Find the peak in omega [0,pi]

required_indices1 = bsxfun(@plus,[-1 0 1]',max_frequency_bins);  %Extract required DFT outputs
required_indices2 = mod(required_indices1-1,2*N)+1; 
required_indices3 = bsxfun(@plus,required_indices2,(0:trial-1)*2*N); 

R = R(required_indices3); 
R = R.*exp(+j*pi/2/N*(N-1)*(required_indices2-1)); clear required_indices1 required_indices2 required_indices3
Rre = real(R); Rim = imag(R); clear R;

num = [-1 0 1]*Rre; denum = [-1 2 -1]*Rre;
ratio_real=num./denum;

num = [-1 0 1]*Rim; denum = [-1 2 -1]*Rim;
ratio_imag=num./denum; clear num denum

delta0_ini=0.25; %Taylor Series Expansion point for 1st iteration
[deltahat_real,deltahat_imag] = inverse_mappings(ratio_real,ratio_imag,N,max_frequency_bins,delta0_ini,maxiter);                                      

ind = find(abs(deltahat_real)>1 | isnan(deltahat_real)); deltahat_real(ind) = 0; % Get rid off obviously wrong estimates (|delta| > 1)
ind = find(abs(deltahat_imag)>1 | isnan(deltahat_imag)); deltahat_imag(ind) = 0; % generated at very low SNR values

alpha = Rre(2,:).^2./(Rre(2,:).^2 + Rim(2,:).^2); %Calculate fusion coefficient
deltahat_fused = sum([alpha; 1-alpha].*[deltahat_real; deltahat_imag],1);
deltahat_fused = deltahat_fused + max_frequency_bins-1;  %with the units of 2N-DFT bins

deltahat_fused = deltahat_fused/2; %Return results with the units of N-point DFT bins

%%
function [deltahat_real,deltahat_imag] = inverse_mappings(ratio_real,ratio_imag,N,max_frequency_bins,delta_ini,maxiter)

delta0_real = delta_ini; delta0_imag = delta_ini; 
for iter=1:maxiter
    [A1,A1dot] = get_coefs2delta_real_2N([-1 0 1],[-1 0 1],N,max_frequency_bins-1,delta0_real);
    [B1,B1dot] = get_coefs2delta_real_2N([-1 0 1],[-1 2 -1],N,max_frequency_bins-1,delta0_real);
    const_real = A1./B1;
    cN_real = (A1dot.*B1 - A1.*B1dot)./B1.^2;
    deltahat_real = (ratio_real-const_real)./cN_real + delta0_real;
    %est_real(iter,:) = deltahat_real;
    delta0_real = deltahat_real;
    
    [A2,A2dot] = get_coefs2delta_imag_2N([-1 0 1],[-1 0 1],N,max_frequency_bins-1,delta0_imag);
    [B2,B2dot] = get_coefs2delta_imag_2N([-1 0 1],[-1 2 -1],N,max_frequency_bins-1,delta0_imag);
    const_imag= A2./B2;
    cN_imag= (A2dot.*B2 - A2.*B2dot)./B2.^2;
    deltahat_imag = (ratio_imag-const_imag)./cN_imag + delta0_imag;
    %est_imag(iter,:) = deltahat_imag;
    delta0_imag = deltahat_imag;
 end;

function [Acoef,Adotcoef] = get_coefs2delta_real_2N(indices,factors,N,kp,delta0_real);

cscv = @(l,delta0,kp) sin(pi/2.*(-l+delta0))./sin(pi/2/N.*(-l+delta0)) + sin(pi/2.*(2*kp+l+delta0))./sin(pi/2/N.*(2*kp+l+delta0));
cscdotv = @(l,delta0,kp) pi/2*cos(pi/2*(-l+delta0))./sin(pi/2/N*(-l+delta0))-pi/2/N*cos(pi/2/N*(-l+delta0)).*sin(pi/2*(-l+delta0))./(sin(pi/2/N*(-l+delta0))).^2+...
    +pi/2*cos(pi/2*(2*kp+l+delta0))./sin(pi/2/N.*(2*kp+l+delta0))-pi/2/N*cos(pi/2/N*(2*kp+l+delta0)).*sin(pi/2*(2*kp+l+delta0))./(sin(pi/2/N*(2*kp+l+delta0))).^2;

Acoef = 0; Adotcoef = 0;
for k=1:length(indices),
    Acoef = Acoef + factors(k).*cscv(indices(k),delta0_real,kp);
    Adotcoef = Adotcoef + factors(k).*cscdotv (indices(k),delta0_real,kp);
end;

%%
function [Acoef,Adotcoef] = get_coefs2delta_imag_2N(indices,factors,N,kp,delta0_imag);

cscv = @(l,delta0,kp) sin(pi/2*(-l+delta0))./sin(pi/2/N*(-l+delta0)) - sin(pi/2*(2*kp+l+delta0))./sin(pi/2/N*(2*kp+l+delta0));
cscdotv = @(l,delta0,kp) pi/2*cos(pi/2*(-l+delta0))./sin(pi/2/N*(-l+delta0))-pi/2/N*cos(pi/2/N*(-l+delta0)).*sin(pi/2*(-l+delta0))./(sin(pi/2/N*(-l+delta0))).^2+...
    -pi/2*cos(pi/2*(2*kp+l+delta0))./sin(pi/2/N*(2*kp+l+delta0))+pi/2/N*cos(pi/2/N*(2*kp+l+delta0)).*sin(pi/2*(2*kp+l+delta0))./(sin(pi/2/N*(2*kp+l+delta0))).^2;

Acoef = 0; Adotcoef = 0;
for k=1:length(indices),
    Acoef = Acoef + factors(k)*cscv(indices(k),delta0_imag,kp);
    Adotcoef = Adotcoef + factors(k)*cscdotv (indices(k),delta0_imag,kp);
end;
