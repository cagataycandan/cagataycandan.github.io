function [a_wFB,b_wFB,a_FB,b_FB] = wFBprediction(in,ARorder)
%function [a_wFB,b_wFB,a_FB,b_FB] = wFBprediction(in,ARorder)
%
% Implements weighted forward-backward prediction method 
% which is the first stage of the paper titled
% "Making Forward-Backward Prediction Perform Like Maximum Likelihood in
%  Gaussian Autoregressive Model Parameter Estimation".
% 
% Inputs: 
%  in : Input vector
%  ARorder : order of the AR system
%
% Outputs: 
%  a_wFB,b_wFB : Outputs of weighted Forward-Backward prediction method
%  a_FB,b_FB: Outputs for the conventional Forward-Backward prediction
%  
%August 2018,
%Cagatay Candan
%

in = in(:); 
N = length(in); 
Nint = N-1;

Af = convm(in(1:Nint),ARorder); Af(1:ARorder-1,:)=[]; Af(end:-1:end-ARorder+2,:)=[];
bf = in(ARorder+1:N);

in = flipud(conj(in));
Ab = convm(in(1:Nint),ARorder); Ab(1:ARorder-1,:)=[]; Ab(end:-1:end-ARorder+2,:)=[];
bb = in(ARorder+1:N);

Amat = [Af; Ab]; 
bvec = [bf; bb]; 
a = inv(Amat'*Amat)*Amat'*bvec; 
a_FB = [1 -a.'];
b_FB = sqrt(xHinvRx(a_FB,in)/N);

dvec = sqrt(1:length(bf));
Amat = [bsxfun(@times,dvec(:),Af); bsxfun(@times,dvec(:),Ab)]; 
bvec = [bsxfun(@times,dvec(:),bf); bsxfun(@times,dvec(:),bb);]; 
a = inv(Amat'*Amat)*Amat'*bvec; 
err = Amat*a - bvec;
a_wFB = [ 1 -a.'];
b_wFB = sqrt(xHinvRx(a_wFB,in)/N);

function out = convm(in,numcol);
out = zeros(length(in)+numcol-1,numcol);

for ind=1:numcol,
    out(ind:ind+length(in)-1,ind) = in; 
end;

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
    
%%%
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


    