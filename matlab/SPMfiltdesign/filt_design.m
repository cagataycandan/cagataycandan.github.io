% This script implements magnitude-only least-squares filter design
% as described in Example 2.
%
% Requires: spectral_fact.m  --- spectral factorization (also provided)
%
% September 2018, 
% Cagatay Candan
%
% 

N=16; 
Lpoints = max(128,N*10);          %Number of points in the frequency grid           
wp1 = 0.1*2*pi; wp2 = 0.35*2*pi;  %passband is the interval of [wp1,wp2]
                                    % The following is assumed for cut-off freqs.: 
                                    %       wp1 < wp2, 
                                    %      -0.5*pi < wp1 < 0.5*pi, 
                                    %      -0.5*pi < wp2 < 0.5*pi, 
dum = 2*pi*((0:Lpoints-1)/Lpoints)'; %frequency grid for optimization
w = [dum(Lpoints/2+1:end)-2*pi; dum(1:Lpoints/2)]; %shifted frequency grid
                                                   % (just to be compatible with
                                                   % fftshift function)

%Set desired ESD characteristic
ind = find((w >= wp1) & (w <= wp2));
Rdesired = zeros(Lpoints,1) + 1e-1;
Rdesired(ind) = 1;    

%%

A = exp(-j*kron(w,[1:N-1]));  %Construct A matrix 
vec1 = ones(Lpoints,1);       %Construct vector of all ones

M = [vec1'*vec1    vec1'*A       vec1'*conj(A); ...
     A'*vec1       A'*A          A'*conj(A); ...
     A.'*vec1      A.'*A         A.'*conj(A)];
Mb = [vec1'; A'; A.']*Rdesired;

r = inv(M'*M)*M'*Mb; %Auto-correlation design is completed 

%Check whether designed auto-correlation is Hermitian symmetric or not
if all(abs( r(2:N) - conj(r(N+1:end)) ) < 1e-10), 
    disp('Indeed, designed auto-correlation sequence is Hermitian symmetric.');
else
    disp('Something went wrong, designed auto-correlation sequence is NOT Hermitian symmetric.');
end;

%CALCULATE ENERGY SPECTRAL DENSITY USING TWO METHODS 
dum = vec1*r(1) + 2*real(A*r(2:N)); %Using the constructed equation system 
dum2 = fft([r(1:N); zeros(Lpoints-2*N+1,1);r(end:-1:N+1)]); %ESD by FFT

%Double check whether dum and dum2 are identical or not
if all(abs(dum-fftshift(dum2))<1e-10),
    disp('ESD calculations by the equation system and FFT are identical.');
else
    disp('Somethign went wrong, ESD calculations by the equation system and FFT are NOT identical.');
end;

figure(1), 
subplot(311),
plot(w/2/pi,Rdesired,'linewidth',3); 
hold on; set(gca,'fontsize',11);
plot(w/2/pi,real(dum),'linewidth',2.5); %Plot energy spectral density
hold off; 

xlabel('normalized frequency = \omega / (2 \pi)');
ylabel('R_d(e^{j\omega}) and R_h(e^{j\omega})');
title('Least Squares Approximation to Desired Magnitude Response'); 
legend('Desired response',['N = ' num2str(N) ' coefficient design'],'location','best');
grid on; set(gca,'xtick',-0.5:0.1:0.5);

if any(real(dum)<0), %CHECK energy spectral density is positive valued or not
    disp('Opps negative valued energy spectrum!'); 
    disp('Can''t execute spectral factorization... Sorry...'); 
    return; 
end;

%%
h = spectral_fact(r(1:N));  %Spectral Factorization
subplot(312), 
stem(0:N-1,real(h),'linewidth',2.5); hold on;
stem(0:N-1,imag(h),'linewidth',2.5); hold off; 
set(gca,'fontsize',11); grid on; 
xlabel('sample index'); ylabel('h[n]');
legend('real part','imag part'); title('Filter Impulse Response'); 

subplot(313),
plot(w/2/pi,fftshift(abs(fft(h,Lpoints)).^2),'linewidth',2.5);
set(gca,'fontsize',11);
xlabel('normalized frequency = \omega / (2 \pi)');
ylabel('|H(e^{j\omega})|^2');
title('Filter Magnitude Response'); 
grid on; set(gca,'xtick',-0.5:0.1:0.5);
