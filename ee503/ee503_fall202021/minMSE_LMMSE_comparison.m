MCnum = 10e3; 

xvec = rand(2,MCnum) - 0.5; %Elements of xvec unif. in [-0.5, 0.5]

%Generate binary valued r.v. taking values -1 and 1 with equal probability
dum = rand(1,MCnum); bin_rv = dum<0.5 + 0; %r.v. taking values {0,1}
bin_rv = 2*bin_rv - 1; %r.v. taking values {-1,1}
xvec = xvec + [0.5 0.5]'*bin_rv; %random vector uniform dist. in 1st and 3rd quadrant

%Show realizations of the random vectors 
figure(1),
plot(xvec(1,:),xvec(2,:),'.','markersize',10);
title([num2str(MCnum) ' realizations of random vectors' ]);

%%
x_component = xvec(1,:); y_component = xvec(2,:);
%min.MSE estimator
y_est_minMSE = 0.5*sign(x_component);
error  = y_component - y_est_minMSE; 
bias_minMSE = mean(error); MSE_minMSE = mean(error.^2);
bias_minMSE_th = 0; MSE_minMSE_th =1/12; 
hist(error,100); title(['min.MSE estimator (MCnum =', num2str(MCnum) ')' char(10) ...
                        'bias-MC = ' num2str(bias_minMSE) ', MSE-MC = ' num2str(MSE_minMSE) char(10) ...
                        'bias-theory = ' num2str(bias_minMSE_th) ', MSE-theory = ' num2str(MSE_minMSE_th)]);
disp(['min.MMSE estimator: 0.5*sgn(x)']);                    
%%                     
%Linear MMSE estimator 
% y = w1*x 
Rx = 1/3; ryx = 1/4; wvec = inv(Rx)*ryx; 
y_est_LMMSE = wvec*x_component;
error  = y_component - y_est_LMMSE; 
bias_LMMSE = mean(error); MSE_LMMSE = mean(error.^2);
bias_LMMSE_th = 0; MSE_LMMSE_th = 1/3 - wvec'*ryx; 
hist(error,100); title(['LMMSE estimator (MCnum =', num2str(MCnum) ')' char(10) ...
                        'bias-MC = ' num2str(bias_LMMSE) ', MSE-MC = ' num2str(MSE_LMMSE) char(10) ...
                        'bias-theory = ' num2str(bias_LMMSE_th) ', MSE-theory = ' num2str(MSE_LMMSE_th)]);
str1 = ['yhat = ' num2str(wvec) 'x '];
disp(['LMMSE estimator: ' str1]);

                  
%%                     
%Linear & Cubic - MMSE estimator 
% y = w1*x + w3*x^3
Rx = [1/3 1/5; 1/5 1/7]; ryx = 0.5*[1/2 1/4]'; wvec = inv(Rx)*ryx; 
y_est_LCMMSE = wvec'*[x_component; x_component.^3];
error  = y_component - y_est_LCMMSE; 
bias_LCMMSE = mean(error); MSE_LCMMSE = mean(error.^2);
bias_LCMMSE_th = 0; MSE_LCMMSE_th = 1/3 - wvec'*ryx; 
hist(error,100); title(['LC-MMSE estimator (MCnum =', num2str(MCnum) ')' char(10) ...
                        'bias-MC = ' num2str(bias_LCMMSE) ', MSE-MC = ' num2str(MSE_LCMMSE) char(10) ...
                        'bias-theory = ' num2str(bias_LCMMSE_th) ', MSE-theory = ' num2str(MSE_LCMMSE_th)]);

str1 = ['yhat = ' num2str(wvec(1)) 'x  ' num2str(wvec(2)) 'x^3'];
disp(['LC-MMSE estimator: ' str1]); 
%%
%Linear & Cubic & Quintic - MMSE estimator 
% y = w1*x + w3*x^3 + w5*x^5
Rx = [1/3 1/5 1/7
      1/5 1/7 1/9
      1/7 1/9 1/11]; 
ryx = 0.5*[1/2 1/4 1/6]'; 
wvec = inv(Rx)*ryx; 
y_est_LCQMMSE = wvec'*[x_component; x_component.^3 ; x_component.^5];
error  = y_component - y_est_LCQMMSE; 
bias_LCQMMSE = mean(error); MSE_LCQMMSE = mean(error.^2);
bias_LCQMMSE_th = 0; MSE_LCQMMSE_th = 1/3 - wvec'*ryx; 
hist(error,100); title(['LCQ-MMSE estimator (MCnum =', num2str(MCnum) ')' char(10) ...
                        'bias-MC = ' num2str(bias_LCQMMSE) ', MSE-MC = ' num2str(MSE_LCQMMSE) char(10) ...
                        'bias-theory = ' num2str(bias_LCQMMSE_th) ', MSE-theory = ' num2str(MSE_LCQMMSE_th)]);

str1 = ['yhat = ' num2str(wvec(1)) 'x  ' num2str(wvec(2)) 'x^3 +' num2str(wvec(3)) 'x^5'];
disp(['LCQ-MMSE estimator: ' str1]);
%%
%%
%Linear & Cubic & Quintic & Septic- MMSE estimator 
% y = w1*x + w3*x^3 + w5*x^5 + w7*x^7
Rx = [1/3 1/5 1/7 1/9
      1/5 1/7 1/9 1/11
      1/7 1/9 1/11 1/13
      1/9 1/11 1/13 1/15]; 
ryx = 0.5*[1/2 1/4 1/6 1/8]'; 
wvec = inv(Rx)*ryx; 
y_est_LCQSMMSE = wvec'*[x_component; x_component.^3 ; x_component.^5 ; x_component.^7];
error  = y_component - y_est_LCQSMMSE; 
bias_LCQSMMSE = mean(error); MSE_LCQSMMSE = mean(error.^2);
bias_LCQSMMSE_th = 0; MSE_LCQSMMSE_th = 1/3 - wvec'*ryx; 
hist(error,100); title(['LCQS-MMSE estimator (MCnum =', num2str(MCnum) ')' char(10) ...
                        'bias-MC = ' num2str(bias_LCQSMMSE) ', MSE-MC = ' num2str(MSE_LCQSMMSE) char(10) ...
                        'bias-theory = ' num2str(bias_LCQSMMSE_th) ', MSE-theory = ' num2str(MSE_LCQSMMSE_th)]);

str1 = ['yhat = ' num2str(wvec(1)) 'x  ' num2str(wvec(2)) 'x^3 +' num2str(wvec(3)) 'x^5 ' num2str(wvec(4)) 'x^7'];
disp(['LCQS-MMSE estimator: ' str1]); 

%%
figure(2), 
dumx = linspace(-1,1,101);
plot(dumx, 0.5*sign(dumx),'--'); hold all;

ord_highest = 13; %Highest polynomial degree of parametric estimator
for this_order = 1:2:ord_highest,
    matdim = (this_order - 1)/2 + 1; 
    first_row = (1:2:this_order) + 2;
    dum_mat = ones(matdim,1)*first_row + (0:matdim-1)'*2*ones(1,matdim);
    Rx = 1./dum_mat; 
    ryx = 0.25./(1:matdim)';
    wopt = inv(Rx)*ryx;
    min_error_poly_est(matdim) = 1/3 - wopt'*ryx;
    
    dum = [wopt'; zeros(1,matdim)]; dum = [0; dum(:)]; 
    polycoef_vec = flipud(dum);
    str1 = poly2str(polycoef_vec,'x');
    disp(['Estimator #' num2str(matdim) ': ' str1]);
    
    figure(2), plot(dumx, polyval(polycoef_vec,dumx));
end;
figure(2), 
dum =  [{'minMSE'} legendstr('Estimator Order:',1:2:ord_highest)];
title('Estimators');
legend(dum,'location','SouthEast'); grid on;
hold off;
%%
figure(3),
plot(1:2:ord_highest,min_error_poly_est,'o-'); hold on;
minMSE = 1/12; plot(1:2:ord_highest,minMSE*ones(size(min_error_poly_est)),'--'); 
hold off;
dum =  {'Polynomial Estimators','minMSE estimator',};
legend(dum,'location','NorthEast'); grid on;
xlabel('Polynomial order'); ylabel('Estimator MSE')
set(gca,'xtick',1:2:ord_highest);