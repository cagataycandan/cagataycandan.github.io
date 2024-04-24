function [anewt,bnewt] = approx_MLc(in,a0)
%function [anewt,bnewt] = approx_MLc(in,a0)
%
% Implements the second stage of the proposed method in the paper titled
% "Making Forward-Backward Prediction Perform Like Maximum Likelihood in
%  Gaussian Autoregressive Model Parameter Estimation".
% 
% Inputs: 
%  in : Input vector
%  a0 : Initial estimate for a vector
%
% Outputs: 
%  [anewt,bnewt] : Outputs take maximize the likelihood around the initial
%  condition
%  
%August 2018,
%Cagatay Candan
%

a0 = a0(:).';
ARorder=size(a0,2)-1;
iternum = 3; %Always 3 iterations!

for iter=1:iternum,
a0 = a0(2:end).';
circA = xHinvRx([1 a0.'],in); 
%%
N = length(in); 
B1T = hankel(conj(in),[conj(in(end)) zeros(1,ARorder)]);
B2H = toeplitz([in(end); zeros(ARorder-1,1)],in(end:-1:end-ARorder+1));

b1 = B1T*[1; a0(:)]; M1 = B1T(:,2:end);
b2 = B2H*a0(:); M2 = B2H; 

gammas0 = atog([1 a0.']); gammas0 = gammas0(:);
[xx,xxx,J,Jc]=jacobian_a_wrt_gamma(gammas0);
J =  J(2:end,:); Jc =  Jc(2:end,:); 

bM21 = J'*(M1'*M1 - M2'*M2)*J + conj(Jc'*(M1'*M1-M2'*M2)*Jc);
bM22 = J'*(M1'*M1 - M2'*M2)*Jc + conj(Jc'*(M1'*M1-M2'*M2)*J);
bb2 = J'*(M1'*b1 - M2'*b2) + conj(Jc'*(M1'*b1 - M2'*b2)); 

bM = 1/circA*[conj(bM22) conj(bM21); bM21 bM22]; 
bb = 1/circA*[conj(bb2); bb2];

%%
fzc = @(gammas) -gammas./(1 - abs(gammas).^2);
fzczc = @(gammas) (-gammas.^2)./(1 - abs(gammas).^2).^2;
fzcz =  @(gammas) -1./(1 - abs(gammas).^2).^2;

dumvec = -1/N*(1:ARorder).';
bv2 = fzc(gammas0).*dumvec;  
bN21 = diag( fzcz(gammas0).*dumvec );
bN22 = diag( fzczc(gammas0).*dumvec ); 

bN = [conj(bN22) conj(bN21); bN21 bN22]; 
bv = [conj(bv2); bv2];

%%
Mt = bM + bN; 
vt = bb + bv; 
if isreal(in)==2,
    delta_gammas = inv(Mt(1:ARorder,1:ARorder)+Mt(1:ARorder,ARorder+1:end))*(-vt(1:ARorder,:));
else 
    delta_gammas = inv(Mt)*(-vt);
    delta_gammas = delta_gammas(1:ARorder); 
end;

anewt = gtoa(gammas0+delta_gammas).'; 
bnewt = sqrt(xHinvRx(anewt.',in)/N);
a0 = anewt; 
end;

%%%%%%%%%
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
%%%%%%%%%
%%%%%%%%%
%%%%%%%%%
function a=gtoa(gamma)

%GTOA	Step-up recursion

%----

%USAGE	a = gtoa(gamma) 

%

%	The direct-form filter coefficients

%		a=[1 a(1) ... a(p)]

%	are derived from the reflection coefficients gamma

%	using the step-up recursion.

%

%  see also ATOG, ATOR, GTOR, RTOA, RTOG

%

%---------------------------------------------------------------

% copyright 1996, by M.H. Hayes.  For use with the book 

% "Statistical Digital Signal Processing and Modeling"

% (John Wiley & Sons, 1996).

%---------------------------------------------------------------



a=1;

gamma=gamma(:);

p=length(gamma);

for j=2:p+1;

	a=[a;0] + gamma(j-1)*[0;conj(flipud(a))];

end
%%%%%%%%%
%%%%%%%%%
%%%%%%%%%
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
%%%%%%%%%
%%%%%%%%%
%%%%%%%%%

function [Jr,Ji,G,Gc]=jacobian_a_wrt_gamma(gamma)
%Returns Jacobian matrix of (a_p wrt to gamma1, gamma2, ... , gammaP
%
% %Example:  
%     gammas = [0.3 0.2j -0.2]; 
%     [Jr,Ji,G,Gc] = jacobian_a_wrt_gamma(gammas);
% 
%     delta = rand(size(gammas))*0.003 +j*rand(size(gammas))*0.0001;
%     a1 = gtoa(gammas+delta), 
%     a2 = gtoa(gammas) + Jr*real(delta(:)) + Ji*imag(delta(:)),
%     a3 = gtoa(gammas) + G*delta(:) + Gc*conj(delta(:)),


a=1;

gamma=gamma(:);

p=length(gamma);
a_vecs = cell(1,p); 

for jind=2:p+1;
	a=[a;0] + gamma(jind-1)*[0;conj(flipud(a))];
    a_vecs{jind-1} = a; 
end

Jr = zeros(p+1,p); Ji = zeros(p+1,p); 

for gammaind=1:p,
    if gammaind==1, vec = 1; else vec = a_vecs{gammaind-1}; end;  
    vecr = [0; conj(flipud(vec))];
    veci = 1i*vecr;
    for jind=(gammaind+1):p; 
     vecr=[vecr;0] + gamma(jind)*[0;conj(flipud(vecr))];
     veci=[veci;0] + gamma(jind)*[0;conj(flipud(veci))];
    end;
Jr(:,gammaind) = vecr;     
Ji(:,gammaind) = veci;
end;
%Ji = Ji; 

G  = Jr/2 - 1i*Ji/2/1;
Gc = Jr/2 + 1i*Ji/2/1;
%%%%