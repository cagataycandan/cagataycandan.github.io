function [min_jnr, Pb_attained] = min_jnr_vs_beta(P_tb, Pb_tolerance, Pb_classic_minimum, varargin)
%% function [min_jnr, Pb_attained] = min_jnr_vs_beta(P_tb, Pb_tolerance, Pb_classic_minimum, snr_db, beta_db, omega_db)
%
% This function finds the minimum required JNR (min_jnr) for classical SLB
% system that achieves a minimum $P_b$ (Pb_attained) within a tolerance with
% respect to optimum SLB detector. It plots min_jnr and Pb_attained with 
% respect to different antenna gain margins $\beta^2$ (beta_db).
%
% Input Parameters:
% omega_db    : Side lobe gain of main antenna, 1x1 scalar, in dB, default 
% value: -30 dB.
% beta_db     : Antenna gain margin between axuilarry antenna and sidelobe of
% main antenna, 1xM vector, in dB, default value: 5:1:10.
% snr_db      : SNR of interest, 1xM vector, in dB, default value: [9 12 15]. 
% P_tb        : Desired probability of target blanking, ex: 0.05, required.
% Pb_tolerance: Tolerance with respect to optimum SLB detector, ex: 0.05.
% Pb_classic_minimum : Desired minimum jammer blanking probability
%
% Sample Run: 
% P_tb, Pb_tolerance and Pb_classic_minimum are required to input to run the
% function, snr_db, beta_db and omega_db have default values as defined
% above.
%
% To generate Fig. 8, run >> min_jnr_vs_beta(0.05, 0.05, 0.90); 
% To generate Fig. 8 with SNR of [10 15 20] dB 
% run >> min_jnr_vs_beta(0.05, 0.05, 0.90, [10 15 20]);
%
% To generate Fig. 8 with P_tb = 0.1, Pb_tolerance = 0.05, 
% Pb_classic_minimum = 0.85, SNR = [10 15 20], $\beta^2$ = [4:2:20] dB and 
% $\omega^2$ = -40 dB,
% run >> min_jnr_vs_beta(0.1, 0.05, 0.85, [10 15 20], [4:2:20], -40);
%
% Caution: SNR should be high enough such that requiered $F$ will be less
% than chosen $\beta^2$. If $F$ is not less than $\beta^2$, then $P_b$ of 0.90 will
% not be attained and infinite loop can occur. To avoid this, first find
% the minimum SNR that achieves $F << \beta^2$.
%
% Osman Coskun
% Nov. 2015 
%
numvarargs = length(varargin);
if numvarargs > 3
    error('find_min_jnr_wrt_beta:TooManyInputs', ...
        'requires at most 3 optional inputs');
end
default_values = {[9 12 15], [5:10], -30};
default_values(1:numvarargs) = varargin;
[snr_db, beta_db, omega_db] = default_values{:};
%initialization
Pb_attained = zeros (length(snr_db), length(beta_db));
Pb_new = ones (length(snr_db), length(beta_db));
min_jnr = Pb_new;
aux_legend = cell (size(snr_db));
for ii = 1: length (snr_db)    
    [F_db,~] = compute_F_SW1(P_tb, snr_db (ii), omega_db);
    for jj = 1:length(beta_db)
        jnr_db = -33;
        while ~((Pb_new (ii,jj) - Pb_attained (ii,jj) <= Pb_tolerance )&& Pb_attained (ii,jj)>=Pb_classic_minimum)           
            Q = compute_Q_matrix (snr_db (ii), jnr_db, omega_db, beta_db(jj));
            [a, b] = compute_hypothesis_a_b ('H1', Q, snr_db (ii), jnr_db, omega_db, beta_db(jj));
            [eta, ~] = compute_threshold (P_tb, a, b);    
            [a, b] = compute_hypothesis_a_b ('H2', Q, snr_db (ii), jnr_db, omega_db, beta_db(jj));
            Pb_new (ii,jj) = compute_Pb_new (a, b, eta);
            Pb_attained (ii,jj) = compute_Pb_classic_SW1 (F_db, jnr_db, beta_db(jj));
            jnr_db = jnr_db + .1;
        end
        min_jnr (ii,jj) = jnr_db;
    end
    aux_legend(ii) = {['SNR = ' num2str(snr_db (ii)) ' dB' ' $F = $ ' sprintf('%1.2f', F_db) ' dB' ]};
end
%% 
figure(2)
h2= plot(beta_db, Pb_attained');
ylabel ('$P_b$ (Maisel)','FontSize', 8,'interpreter', 'latex')
xlabel ('$\beta^2$  (dB)','FontSize', 8,'interpreter', 'latex')
set (gca(), 'xtick', beta_db, 'FontSize', 8);
title (['$\omega^2=$ ' num2str(omega_db) ' dB, ' '$P_{tb}=$ ' num2str(P_tb) ', min.' '$P_{b}=$ ' num2str(Pb_classic_minimum)],'interpreter', 'latex');
grid on
legend(aux_legend, 'Location','northeast','FontSize', 8,'interpreter', 'latex');
figure(1)
h1 = plot(beta_db, min_jnr');
set (gca(), 'xtick', beta_db, 'FontSize', 8);
grid on
ylabel ('Min. JNR','interpreter', 'latex', 'FontSize', 8)
xlabel ('$\beta^2$  (dB)','interpreter', 'latex','FontSize', 8)
title (['$\omega^2=$ ' num2str(omega_db) ' dB, ' '$P_{tb}=$ ' num2str(P_tb) ', min.' '$P_{b}=$ ' num2str(Pb_classic_minimum)],'interpreter', 'latex');
marker_style =  ['s', '<', 'o', '+', '>', 'v', ];
marker_face_color = ['r', 'b', 'g', 'm'];
for jj = 1: length (snr_db) 
    set(h1(jj),'Linestyle', '-', 'Marker', marker_style (jj), 'Linewidth', 0.6, ...
        'MarkerSize',4, 'MarkerFaceColor','w', 'MarkerEdgeColor','k','Color',marker_face_color(jj)); 
    set(h2(jj),'Linestyle', '-', 'Marker', marker_style (jj), 'Linewidth', 0.6, ...
        'MarkerSize',4, 'MarkerFaceColor',marker_face_color(jj), 'MarkerEdgeColor','k','Color',marker_face_color(jj)); 
end
legend(aux_legend, 'Location','northeast','FontSize', 8,'interpreter', 'latex');
% *************
% End of main function
%%
function [F_db,Pb_found] = compute_F_SW1(Pb_desired, jnr_db, beta_db)
beta = db2pow(beta_db);      % beta in linear scale
jnr = db2pow (jnr_db);
tolerance= Pb_desired/1e5;    % Tolerance value for computation of 
initial = 0.000001;            % Initial value for bisection search
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
