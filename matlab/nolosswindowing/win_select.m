function win_index = win_select2(r,wins,alpha_vec)
% r : input vector (column vector)
% wins : matrix of window function, each column contains a window function
%        wins = [win1(:) win2(:) win3(:) ... winN(:)]; 
%        1st column of wins should have the lowest peak-side-lobe level 
%                             (least interference suppressing window)
%        2nd column of wins should have a better PSL than 1st column 
%        ...
%        last column of wins should have the best PSL
%                             (most interference suppressing window)
%
% alpha_vec : DFT bins of interest
% 
%Cagatay Candan
%Oct.2016
%
% thisx=27; dum1 = mywin_select(r(:,thisx),wins); [dum1(7) selected_win2(thisx)]

[N,numwin] = size(wins); 
if size(r,1)~=N, 
    disp('Dimensions of the argument "r" and "wins" do not match'); return; 
end;

if exist('alpha_vec')==0, alpha_vec=0:N-1; end; 

%Step 0
bin_away = get_theta1(wins);
dec_boundaries = find_thresholdsfromROC(wins,bin_away);

%Step 1
for ind=2:numwin,
    M(:,:,ind) = generate_Rjn(N,0,bin_away(ind),0);
end;

win_index = zeros(length(alpha_vec),size(r,2));
for ind=1:size(r,2),
    alphaind=0;
    for thisalpha=alpha_vec,
        alphaind=alphaind+1;
        rdum = exp(-j*2*pi/N*thisalpha*(0:N-1)').*r(:,ind);
        win_index(alphaind,ind) = Steps2to4(rdum,N,numwin,M,dec_boundaries);
    end;
end;
    
%
function win_index = Steps2to4(r,N,numwin,M,dec_boundaries)
%Step 2
for ind=2:numwin,
    d_metric(ind) = real(r'*M(:,:,ind)*r)/N;
end;

%Step 3
last_win = numwin;
for ind=2:numwin-1,
    if d_metric(ind) > 2*d_metric(ind+1), last_win = ind; break; end;
end;

d_metric = 10*log10(d_metric);% now in dB

%Step 4
dec_boundaries_in_use = dec_boundaries([1:last_win end]);
dum = d_metric(2);
for win_index=1:(length(dec_boundaries_in_use)-1),
    low_thresh = dec_boundaries_in_use(win_index);
    high_thresh = dec_boundaries_in_use(win_index+1);
    if (dum>low_thresh) & (dum<high_thresh), break; end;
end;
%end of Steps2to4


function bin_away = get_theta1(wins)
%Returns start of stopband for each column of wins matrix
%
[N,numwin] = size(wins);
dum = sum(wins); wins = wins.*repmat(1./dum,N,1);
winsf_m = abs(fft(wins,1024,1));
freq_axis = linspace(0,N-1/1024,1024);

for ind=1:numwin-1,
    dum=find(winsf_m(:,ind+1)<winsf_m(:,ind));
    bin_away(ind) = freq_axis(dum(2));
end;
bin_away = [0 bin_away];
%%% end of function bin_away = get_theta1(wins)

function winselect2_thresholds = find_thresholdsfromROC(wins,bin_away_vec)
% Assumes that windows are ordered such that the least sidelobe supressing 
% window is in the first column of wins, ... , 
% the most sidelobe supressing window is in the last column of wins. 
%

numwin = size(wins,2);
winselect2_thresholds(1)= -500;
for indx=2:numwin,
        dum = get_pd_crossing(wins(:,indx-1),wins(:,indx),bin_away_vec(indx));
        winselect2_thresholds(indx) = dum;  %Select threshold 1 dB less the crossing value
end; %for
winselect2_thresholds(end+1) = 500;
%%%end of function find_thresholdsfromROC

function crossingJNR = get_pd_crossing(win1,win2,bin_away)
N  = length(win1); 
look_bin = 0; %can be any bin number!
s_look = ones(size(win1));  %steering vector for look_bin = 0

Rn = generate_Rjn(N,look_bin,bin_away,0); 
win1 = win1.*s_look;
win2 = win2.*s_look;

A1=abs(win1'*s_look)^2; 
A2=abs(win2'*s_look)^2; 
B1 = sum(abs(win1).^2); 
B2 = sum(abs(win2).^2); 
C1 = real(win1'*Rn*win1);
C2 = real(win2'*Rn*win2);

dum = (A1*B2 - A2*B1)/(A2*C1-A1*C2);
crossingJNR = 10*log10(dum); 
%%%end of function get_pd_crossing


function R = generate_Rjn(N,sbin,bin_away,JNR_dB)
%returns Rjn matrix 
JNR = 10^(JNR_dB/10);
theta1  = (sbin+bin_away)*2*pi/N;
theta2  = (sbin+N-bin_away)*2*pi/N;
dum = [(theta2-theta1) (exp(-j*theta2*(1:N-1))-exp(-j*theta1*(1:N-1)))./(1:N-1)/(-j)]; 
R = JNR*toeplitz(dum)/2/pi;







