N=16;            %Number of input samples
SNR_dB = [5 8]; %SNR Values
bins = [4.3 8.25]; %Signal frequency in bins

%Definitions
svec = @(bin) exp(j*2*pi/N*bin*(0:N-1)');
SNR = 10.^(SNR_dB/10);

%windows
wins(:,1) = window(@rectwin,N)'; 
wins(:,2) = window(@hamming,N)'; 
wins(:,3) = chebwin(N,120)';
dum = sum(wins);
wins = wins.*repmat(1./dum,N,1); %now sum(wins) = 1 for all windows
win_string = 'RHC';

%Construct Input for Monte Carlo Run
r = sqrt(SNR(1))*svec(bins(1)) + sqrt(SNR(2))*svec(bins(2)) + ... 
       1/sqrt(2)*(randn(N,1) + j*randn(N,1));
   
   
%Window Selection
DFTbins = 5:N-1;
out = win_select(r,wins,DFTbins);

disp('*****')
disp(['Scenario:' char(10) 'N = ' num2str(N) ',' char(10) ...
     '2 signal components with SNR = ' num2str(SNR_dB(1)) ' and ' num2str(SNR_dB(2)) ' dB,' ...
     char(10) ' located at DFT bins ' num2str(bins(1)) ' and ' num2str(bins(2)) '.']);
str1 = sprintf('DFT Bins: \t');
str2 = sprintf('%d\t',DFTbins);
str3 = sprintf('Selected Windows: \t');
str4 = sprintf('%c\t',win_string(out));
str5 = sprintf('%s\n%s\n%s\n%s',str1,str2,str3,str4);
disp(str5);
disp('*****')

