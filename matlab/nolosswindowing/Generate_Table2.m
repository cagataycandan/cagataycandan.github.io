N=16;         %Number of input samples
MCnum = 10e3; %Monte Carlo Run number 
SNR_dB = [0 35]; %SNR Values
bins = [0.1 6.25]; %Signal frequency in bins

%Definitions
svec = @(bin) exp(j*2*pi/N*bin*(0:N-1)');
SNR = 10.^(SNR_dB/10);

%windows
wins(:,1) = window(@rectwin,N)'; 
wins(:,2) = window(@hamming,N)'; 
wins(:,3) = chebwin(N,120)';
dum = sum(wins);
wins = wins.*repmat(1./dum,N,1); %now sum(wins) = 1 for all windows

%Construct Input for Monte Carlo Run
fading1 = 1/sqrt(2)*(randn(1,MCnum) + j*randn(1,MCnum));
fading2 = 1/sqrt(2)*(randn(1,MCnum) + j*randn(1,MCnum));
Rmat = sqrt(SNR(1))*repmat(svec(bins(1)),1,MCnum).*repmat(fading1,N,1) + ...
       sqrt(SNR(2))*repmat(svec(bins(2)),1,MCnum).*repmat(fading2,N,1) + ...
       1/sqrt(2)*(randn(N,MCnum) + j*randn(N,MCnum));
   
   
%Window Selection
DFTbins = 0:N-1;
out = win_select(Rmat,wins,DFTbins);

%Display results
out1 =  ( sum(out==1,2)/MCnum*100 + eps) ;
out2 =  ( sum(out==2,2)/MCnum*100 - eps) ;
out3 =  ( sum(out==3,2)/MCnum*100 ) ;
[DFTbins; out1'; out2'; out3']