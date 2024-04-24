function Pb = Pb_vs_beta(P_tb, varargin)
%% function Pb = Pb_vs_beta(P_tb, jnr_db, snr_db, beta_db, omega_db)
%
% This function finds and plots the probability of blanking jammer (Pb) 
% for optimum and Maisel SLB systems with respect to different antenna 
% gain margins $\beta^2$ (beta_db).
%
% Input Parameters:
% P_tb        : Desired probability of target blanking, ex: 0.05.
% snr_db      : SNR of interest, 1xM vector, in dB, default = [9 10 12] dB.
% jnr_db      : JNR of interesr, scalar, in dB, default = 5 dB.
% beta_db     : Antenna gain margin between axuilarry antenna and sidelobe of
% main antenna, 1xM vector, in dB, default = [5:1:10] dB.
% omega_db    : Side lobe gain of main antenna, 1x1 scalar, default = -30 dB.
%
% Sample Run: 
% P_tb is required to input to run the function. jnr_db, snr_db,  beta_db and omega_db 
% have default values as defined above.
%
% To generate Fig. 7, run >> Pb_vs_beta(0.05); 
%
% To generate Fig. 7 with P_tb = 0.1, JNR = 8 dB, SNR = [5 10 15],  
% $\beta^2$ = [2:2:20] dB and $\omega^2$ = -40 dB,
% run >> Pb_vs_beta (0.1, 8, [5 10 15], [2:2:20], -40);
%
% Osman Coskun
% Nov. 2015
%
numvarargs = length(varargin);
if numvarargs > 4
    error('plot_Pb_wrt_beta:TooManyInputs', ...
        'requires at most 3 optional inputs');
end
default_values = {5, [9 12 15], [5:10], -30};
default_values(1:numvarargs) = varargin;
[jnr_db, snr_db, beta_db, omega_db] = default_values{:};

main_legend = cell (2*length(snr_db),1);
Pb = zeros (2*length(snr_db),length(beta_db));
for ii = 1: length (snr_db)    
    [F_db,~] = compute_F_SW1(P_tb, snr_db (ii), omega_db);
    for jj = 1:length(beta_db)
        Q = compute_Q_matrix (snr_db (ii), jnr_db, omega_db, beta_db(jj));
        [a, b] = compute_hypothesis_a_b ('H1', Q, snr_db (ii), jnr_db, omega_db, beta_db(jj));
        [eta, ~] = compute_threshold (P_tb, a, b);    
        [a, b] = compute_hypothesis_a_b ('H2', Q, snr_db (ii), jnr_db, omega_db, beta_db(jj));
        Pb (2*ii-1,jj) = compute_Pb_new (a, b, eta);
        Pb (2*ii,jj) = compute_Pb_classic_SW1 (F_db, jnr_db, beta_db(jj));
    end
    main_legend  (2*ii-1) = {['Optimum, SNR = ' num2str(snr_db (ii)) ' dB']};
    main_legend (2*ii) = {['Maisel, SNR = ' num2str(snr_db (ii)) ' dB' ' $F =  $ ' sprintf('%1.2f', F_db) ' dB']};
end
% figure(1)
%%
h1 = plot(beta_db, Pb');
y_tick = 0:0.05:1;
set (gca(), 'xtick', beta_db, 'FontSize', 8);
set(gca(), 'ytick', y_tick, 'FontSize', 8);
grid on
ylabel ('$P_b$','Interpreter','latex')
xlabel ('$\beta^2$  (dB)', 'Interpreter','latex')
title ([ '$\omega^2=$ ' num2str(omega_db) ' dB, ' 'JNR = ' num2str(jnr_db)... 
    ' dB, ' '$P_{tb}=$ ' num2str(P_tb)],'FontSize', 8, 'interpreter', 'latex');
marker_style =  ['s', '<', 'o', '+', '>', 'v', ];
for jj = 1: length (snr_db)
    set(h1(2*jj-1),'Linestyle', '-' , 'Marker', marker_style (jj), 'Linewidth', 0.5, ...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','w',...
        'MarkerSize',3, ...
        'Color', 'Red');       
    set(h1(2*jj),'Linestyle', ':', 'Marker', marker_style (jj), 'Linewidth', 0.5, ... 
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','w',...
        'MarkerSize',3, ...
        'Color', 'Blue');   
end
legend(main_legend, 'Location','southeast','FontSize', 8,'interpreter', 'latex');
set(gca, 'TickLabelInterpreter','latex');
% End of main function
%%
function [F_db,Pb_found] = compute_F_SW1(Pb_desired, jnr_db, beta_db)
beta = db2pow(beta_db);      % beta in linear scale
jnr = db2pow (jnr_db);
tolerance= Pb_desired/1e5;    % Tolerance value for computation of 
initial = 0.000001;            % ?nitial value for bisection search
final= 200;               % Final value for bisection search
Pb = @(F) 0.5.*(1- (1./(1+F)).*(jnr.*(F-beta)-F-1)./sqrt((jnr.*(F-beta)+F+1).^2+ 4.*jnr.*(F+1).*beta) - (F./(1+F)).*(jnr*(F-beta)+F+1)./sqrt((jnr*(F-beta)-F-1).^2+ 4*jnr*(F+1).*F));
initial_Pb = Pb(initial);
error = initial_Pb - Pb_desired;
while abs(error) > tolerance;
    F_linear = 0.5*(initial + final);
    if Pb(F_linear) >= Pb_desired
        initial = F_linear;
    else
        final = F_linear;
    end
    error = Pb(F_linear) - Pb_desired;
end
F_db = 10*log10(F_linear);
Pb_found = Pb(F_linear);
%*****************
function [Q] = compute_Q_matrix (snr_db, jnr_db, omega_db, beta_db)
%[Q] = compute_Q_matrix (snr_db, jnr_db, omega_db, beta_db)
beta = db2pow(beta_db);        % beta in linear scale
omega = db2pow(omega_db);      % omega in linear scale 
jnr = db2pow (jnr_db);         % JNR in  linear scale
snr = db2pow (snr_db);         % SNR in linear scale
% Compute covariance matrices
C1 = [ (snr +1), sqrt(omega)*snr; sqrt(omega)*snr, (omega*snr +1)]; 
C2 = [ (jnr +1), sqrt(beta)*jnr; sqrt(beta)*jnr, (beta*jnr +1)];
Q = inv(C1) - inv(C2);
%*****************
function [a, b] = compute_hypothesis_a_b (hypothesis, Q, snr_db, jnr_db, omega_db, beta_db)
% [a, b] = compute_hypothesis_a_b (hypothesis, Q, snr_db, jnr_db, omega_db, beta_db)
beta = db2pow(beta_db);        % beta in linear scale
omega = db2pow(omega_db);      % omega in linear scale 
jnr = db2pow (jnr_db);         % JNR in  linear scale
snr = db2pow (snr_db);         % SNR in linear scale
% Compute covariance matrices
C1 = [ (snr +1), sqrt(omega)*snr; sqrt(omega)*snr, (omega*snr +1)]; 
C2 = [ (jnr +1), sqrt(beta)*jnr; sqrt(beta)*jnr, (beta*jnr +1)];
if hypothesis == 'H1'
    R = C1/2;
else
    R = C2/2;
end
r = sum(sum(R.*Q))/(-4*det(R)*det(Q));
a = sqrt(r^2+1/(-4*det(R)*det(Q))) - r ;
b = sqrt(r^2+1/(-4*det(R)*det(Q))) + r ;
%*****************
function [eta, threshold] = compute_threshold (P_b, a, b)
threshold = b/(a+b);
if P_b <= threshold
    eta = -(log((P_b*(a + b))/b) )/a;
else
    eta = (log(-((a + b)*(P_b - 1))/a) )/b;
end
%*****************
function Pb_temp = compute_Pb_new (a_temp, b, eta_temp)
f1 = @(x)(b/(a_temp+b)*exp(-a_temp.*x));
f2 = @(x)(a_temp/(a_temp+b))*(1-exp(b.*x))+b/(a_temp+b);
Pb_temp = zeros (size(eta_temp));
for i = 1: length (eta_temp)
    if eta_temp(i) >= 0
        Pb_temp(i) = f1 (eta_temp(i));
    else
        Pb_temp(i) = f2 (eta_temp(i));
    end
end
%*****************
function Pb = compute_Pb_classic_SW1 (F_db, jnr_db, beta_db)
beta = db2pow(beta_db);        % beta in linear scale
jnr = db2pow (jnr_db);         % JNR in  linear scale
F = db2pow (F_db);             % F in linear scale
Pb_formula = @(F, jnr, beta) 0.5.*(1- (1./(1+F)).*(jnr.*(F-beta)-F-1)./sqrt((jnr.*(F-beta)+F+1).^2+ 4.*jnr.*(F+1).*beta) - ...
(F./(1+F)).*(jnr*(F-beta)+F+1)./sqrt((jnr*(F-beta)-F-1).^2+ 4*jnr*(F+1).*F));
Pb = Pb_formula (F, jnr, beta);
