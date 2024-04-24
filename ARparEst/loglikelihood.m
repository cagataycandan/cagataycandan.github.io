function out = loglikelihood(a,b,x)
x = x(:); 
N = length(x); 
%if N>400, disp('x should have a smaller dimension'); return; end;
%r = my_ator(a,1,N);
%Rf = toeplitz(r,r');
sigmasq = abs(b)^2; 
dum1 = xHinvRx(a,x);

gammas = atog(a);
dum2 = 0.5*sum(log((1-abs(gammas).^2)).*(1:length(gammas)));

%out = -0.5*log(det(Rf)) - N/2*log(sigmasq) -1/2/sigmasq*real(x'*inv(Rf)*x); 
out = dum2 - N/2*log(sigmasq) -1/2/sigmasq*real(dum1); 


function output = xHinvRx(a,x)
%function output = xHinvRx(a,x)
% Generates quadratic form of x'*inv(Rf)*x
% Inputs: 
%       a : [1 a1 a2 ... aP]
%       x : (N+1) dimensional vector 
% Output: 
%       output : x'*inv(Rf)*x.   
%
%Here, 
%
% Rf = [r(0) r(-1)  r(-2)  ... r(-N);  
%       r(1) r(0)   r(-1)  ... r(-N+1)
%       .    .       .         .
%       r(N) r(N-1) r(N-2) ... r(0)]
%
%and r(k) is the auto-covariance values of AR(P) process 
%generated with the application of unit variance white noise 
%at the input of the following filter,
%
%  H(z) = 1/(1 + a1z^{-1} + a2z^{-2} + ... aPz^{-P}.
%
% March 1, 2018
% Cagatay Candan
%
% Note: In some references, such as Stoica's textbook (Spectral Analysis of
% Signals), Rf matrix is defined as the Hermitian of the matrix given
% above. To accomodate Stoica's definition, use 
%         output = xHinvRx(conj(a),x))
% to calculate the quadratic form. 
% 
%
% To test the code: 
% a = [1 0.5j]; x = randn(4,1)+j*1*randn(4,1); 
% r = my_ator(a,1,length(x));Rf = toeplitz(r,r'); 
% [x'*inv(Rf)*x xHinvRx(a,x)]
%
%

gamma = atog(a); 
ARorder = length(a) - 1; 
N = length(x);
sigmasq = 1/N*real(x'*x);
e = x; f = x;
alpha = e(1); beta = f(end);
etilde = e(2:end); ftilde = f(1:end-1);
c = etilde'*ftilde; d = 2*N*sigmasq - abs(alpha)^2 - abs(beta)^2;
hsq = 0; 
for order = 1:ARorder,
    k = gamma(order);
    V = (1 + abs(k)^2)*d + 4*real(k*c);
    hsq = (1 - abs(k)^2)*(hsq + abs(alpha)^2 + abs(beta)^2);
    e = etilde + k*ftilde;
    f = conj(k)*etilde + ftilde;
    alpha = e(1); beta = f(end);
    etilde = e(2:end); ftilde = f(1:end-1);
    c = etilde'*ftilde; 
    d = V - abs(alpha)^2 - abs(beta)^2;
end;
output = (hsq + V)/2;


function gamma=atog(a)

%ATOG	Step-down recursion

%----

%USAGE: gamma=atog(a)

%

%	Finds the reflection coefficients gamma from the

%	direct-form filter coefficients a(k).

%

%  see also ATOR, GTOA, GTOR, RTOA, RTOG

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



a=a(:);

p=length(a);

a=a(2:p)/a(1);

gamma(p-1)=a(p-1);

for j=p-1:-1:2;

	a=(a(1:j-1) - gamma(j)*flipud(conj(a(1:j-1))))./ ...
	  (1 - abs(gamma(j))^2);

	gamma(j-1)=a(j-1);
end

