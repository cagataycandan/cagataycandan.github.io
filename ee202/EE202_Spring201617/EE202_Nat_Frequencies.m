A=[-1.2500 0.0625; 1.0000 -1.2500], % A matrix --> nat. freq = {-1, -1.5}
%A=[-1.5 0.125; 2 -1.5], % A matrix --> nat. freq = {-1, -2}
%A=[-3  0.5; 8 -3], % A matrix --> nat. freq = {-1, -5}
%A=[-5.5  1.125; 18 -5.5], % A matrix --> nat. freq = {-1, -10}

x0 = [1; 2];    %initial condition 
tsim = 10;      %simulation time (sec.) 

[Evec,Eval] = eig(A),
if any(diag(Eval)>=0), disp 'Warning : System is unstable'; end; 

%Solve state equation
%Numerical solution
xdot = @(t,x) A*x;
[t,xnumerical] = ode45(xdot,[0;tsim],x0);

%Analytical solution
coefs = inv(Evec)*x0; 
xanalytical = coefs(1)*Evec(:,1)*exp(Eval(1,1)*t')  + coefs(2)*Evec(:,2)*exp(Eval(2,2)*t');

%Show state variables as a function of time 
figure(1), 
plot(t,xnumerical,'x-'); hold on; 
plot(t,xanalytical,'--','linewidth',3); hold off; 
xlabel('t (sec.)'); 
legend('x_1(t) (numerical)','x_2(t) (numerical)','x_1(t) (analytical)','x_2(t) (analytical)');
set(gca,'fontsize',12),
grid on; 

%Show state trajectory
figure(2), plot(xnumerical(:,1),xnumerical(:,2));  
xlabel('x_1(t)');  ylabel('x_2(t)');
hold on; 

dum = coefs(1)*Evec(:,1); plot([0 dum(1)],[0 dum(2)],'--'); 
dum = coefs(2)*Evec(:,2); plot([0 dum(1)],[0 dum(2)],'--'); 
hold off
h = legend('state trajectory','slow mode eigenvec','fast mode eigenvec','location','NorthWest'); 
title(['Natural Frequencies = [' num2str(Eval(1,1)) ',' num2str(Eval(2,2)) ']']); 
set(gca,'fontsize',12),
grid on;


