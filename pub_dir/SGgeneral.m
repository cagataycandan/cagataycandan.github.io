function [h,SGfilter] = SGgeneral(N,L,b);
% function [h,SGfilter] = SGfracdelay(N,L,delay);
%
% Design Savitzky - Golay "general" filter 
%  with (2N+1) points and L'th order
% 
%
% Input : 
%  N, L : Projects (2N+1) dimensional input to 
%         the L'th order polynomial space (L+1 dim. subspace)
% 
%  b    : (2N+1)x1 dimensional vector whose elements are
%          b(k+1) = desired mapping from polynomial space 
%                   { t^k, k={0,1,..., 2N} } to real line
% 
%  For 
%    smoothing (evaluation at t=0) : b = [1 0 0 0 0 0 ... 0 ]';
%    1st derivative : b = [0 1 0 0 0 0 ... 0]';
%    2nd derivative  : b = [0 0 2 0 0 0 ... 0]'; 
%    fractional advance of 1/2 samples : b = [ 1 1/2 1/4 1/8 ...]';
%    fractional delay of d units (|d|<1): b = [1 (d) (d)^2 (d)^3 ...]';
%    integration : b = [ 1 0 2/3*(1/2)^3 0 2/5*(1/2)^5 0 2/7*(1/2)^7 ...]';
%      
%
% Output : 
%    h : Impulse response of the filter
%    SGfilter : Expression of the filter in terms of finite differences
%
%
%
% August 2013, 
% Cagatay Candan
%

if L>(2*N), disp('Error : L <= 2*N should be satisfied'); return; end;
b = b(:); 

%SPECIAL CASE OF SMOOTHING
if isnumeric(b(1))==1 & b(1)==1 & norm(b(2:end))==0, 
    [h,SGfilter] = SGsmoothingfilter(N,L);
    if nargout==0, 
    disp('SG Filter impulse response (to be used with conv(x,h) command):');
    fliplr(h),
    disp('SG Filter expansion:');
    disp('      D^k ==> \nabla_k defined in the paper');
    disp('      S ==> convolution with all ones, i.e. S==>ones(1,2*N+1)');
    SGfilter,
end;
    return; 
end;
%%END OF SPECIAL CASE

AT = repmat((-N:N),2*N+1,1).^((0:2*N)'*ones(1,2*N+1));
AT = sym(AT);
p = inv(AT'*AT)*AT'*b;  %this is a particular solution!

syms D S; 
dimnull = 2*N-L;

vec=1; Nullmat = zeros(2*N+1,2*N);
Nullmat(N,1)=1/2;Nullmat(N+2,1)=-1/2;
for dum=2:2:2*N,
   vec = conv([1 -2 1],vec); 
   coln = Nullmat(:,dum);
   coln(N+1-dum/2:N+1+dum/2) = vec;
   Nullmat(:,dum) = coln;
   if dum<2*N,
       coln = Nullmat(:,dum+1);
       coln(N-dum/2:N+2+dum/2) = conv(vec,[-1/2 0 1/2]);
       Nullmat(:,dum+1) = coln;
   end;
end;
dum =[ ones(2*N+1,1) Nullmat ]; 
pcoef = inv(dum)*p; 

Nullmat = Nullmat(:,(L+1):end);
Nullmat = sym(Nullmat);

if length(Nullmat)>0,    
    coefs = inv(Nullmat'*Nullmat)*Nullmat'*p;
    dum = pcoef(1)*S + (D*ones(1,2*N)).^[1:(2*N)]*(pcoef(2:end));
    SGfilter = dum - (D*ones(1,dimnull)).^[(L+1):(2*N)]*(coefs);
    %SGfilter = simplify(subs(SGfilter,'D','-D')); %Original 
    SGfilter = simplify(subs(SGfilter,'D',-D)); %Matlab R2022b 
    h = p - Nullmat*coefs; h= h';
else;
    SGfilter = pcoef(1)*S + (D*ones(1,2*N)).^[1:(2*N)]*(pcoef(2:end));
    %SGfilter = simplify(subs(SGfilter,'D','-D')); %Original 
    SGfilter = simplify(subs(SGfilter,'D',-D)); %Matlab R2022b 
    h = p; h= h';
end;

if nargout==0, 
    disp('SG Filter impulse response (to be used with conv(x,h) command):');
    fliplr(h),
    disp('SG Filter expansion:');
    disp('      D^k ==> \nabla_k defined in the paper');
    disp('      S ==> convolution with all ones, i.e. S==>ones(1,2*N+1)');
    SGfilter,
end;


%%%%%%%%%%%%%
function [h,SGfilter] = SGsmoothingfilter(N,L);
% function [h,SGfilter] = SGsmoothingfilter(N,L);
%
% Design Savitzky - Golay smoothing filter 
%  with (2N+1) points and L'th order
%
% (Projects (2N+1) dimensional input to 
%  the L'th order polynomial space)
% 
% Output : 
%    h : Impulse response of the filter
%    SGfilter : Expression of the filter in terms of finite differences
%
% August 2013, 
% Cagatay Candan
%

if L>(2*N), disp('Error : L <= 2*N should be satisfied'); return; end;

p = zeros(2*N+1,1); p(N+1)=1;  % particular solution for smoothing app.
syms D; 
dimnull = 2*N-L;

vec=1; Nullmat = zeros(2*N+1,2*N);
Nullmat(N,1)=1/2;Nullmat(N+2,1)=-1/2;
for dum=2:2:2*N,
   vec = conv([1 -2 1],vec); 
   coln = Nullmat(:,dum);
   coln(N+1-dum/2:N+1+dum/2) = vec;
   Nullmat(:,dum) = coln;
   if dum<2*N,
       coln = Nullmat(:,dum+1);
       coln(N-dum/2:N+2+dum/2) = conv(vec,[-1/2 0 1/2]);
       Nullmat(:,dum+1) = coln;
   end;
end;

Nullmat = Nullmat(:,(L+1):end);
Nullmat = sym(Nullmat);

if length(Nullmat)>0,    
    coefs = inv(Nullmat'*Nullmat)*Nullmat'*p;
    SGfilter = 1 - (D*ones(1,dimnull)).^[(L+1):(2*N)]*coefs;
    h = p - Nullmat*coefs; h= h';
else;
    SGfilter = 1 + D -D;
    h = 1; 
end;

if nargout==0, 
    disp('SG Filter impulse response (to be used with conv(x,h) command):');
    fliplr(h),
    disp('SG Filter expansion:');
    disp('      D^k ==> \nabla_k defined in the paper');
    disp('      S ==> convolution with all ones, i.e. S==>ones(1,2*N+1)');
    SGfilter,
end;

    