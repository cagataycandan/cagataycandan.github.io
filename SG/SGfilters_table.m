N=3;

%Derivative 
fprintf('Derivative SG filters\n');
b = zeros(1,2*N+1); b(2)=1;
for Lthis=1:2*N,
    [h,SGthis]=SGgeneral(N,Lthis,b);
    fprintf('N=%d, L=%d, SG filter =  %s \n',N,Lthis, char(SGthis));
end;
fprintf('----------------------\n');


% Smoothing 
fprintf('Smoothing SG filters\n');
b = zeros(1,2*N+1); b(1)=1;
for Lthis=0:2*N,
    [h,SGthis]=SGgeneral(N,Lthis,b);
    fprintf('N=%d, L=%d, SG filter =  %s \n',N,Lthis, char(SGthis));
end;
fprintf('----------------------\n');

% Fractional Delay
fprintf('Fractional Delay SG filters\n');
syms d; 
b = d.^(0:2*N); b=b(:);
for Lthis=1:2*N,
    [h,SGthis]=SGgeneral(N,Lthis,b);
    fprintf('N=%d, L=%d, SG filter =  %s \n',N,Lthis, char(SGthis));
end;
fprintf('----------------------\n');