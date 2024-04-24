N=16; 
kp_2N = 7; delta_2N = -0.2; %with the units of 2N point DFT bins
trial=10e3; maxiter = 6;
SNR_dB=-10:2:20;

SNR=10.^(SNR_dB./10);
vec = 2*pi/2/N*(kp_2N+delta_2N)*(0:N-1)';
rms_est = zeros(1,length(SNR_dB));
for indSNR=1:length(SNR_dB)
    psi  = 2*pi*rand(1,trial);
    rmat = sqrt(2*SNR(indSNR))*cos(repmat(vec,1,trial)+repmat(psi,N,1))+randn(N,trial);
   
    est_fused = freq_est_cosine_2N(rmat); %est_fused with the units of N point DFT bins
    error = 2*est_fused - (kp_2N+delta_2N);  
    rms_est(:, indSNR) = sqrt(mean(error.^2,2));
end

CRB=12*N/(N*N-1)/((2*pi)^2)./SNR * 4 ; %with the units of 2N point DFT bins
%%

%%
figure(1), clf
semilogy(SNR_dB,rms_est(end,:),'-o','linewidth',2);
set(gca,'fontsize',12);
hold on;
semilogy(SNR_dB,sqrt(CRB),'linewidth',2);
hold off;
legend ('Proposed Method','ACRB or HCRB')
xlabel('SNR(dB)'); ylabel('RMSE (in DFT bins)');
title(['N = ' num2str(N) ', k_p = ' num2str(kp_2N) ', \delta = ' num2str(delta_2N)]);
grid on;

