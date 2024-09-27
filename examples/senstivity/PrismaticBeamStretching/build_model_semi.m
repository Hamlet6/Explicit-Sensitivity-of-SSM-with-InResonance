function [mass,stiff,fnl,DM,DK,dfnl2,dfnl3] = build_model_semi(b,h,L,rho,E,P)

n = 10;
A = b*h;
I = b*h^3/12;
dIdb = h^3/12;
dIdh = b*h^2/4;
mass = rho*A*L*eye(n);
% damp = c*L*eye(n);
tmp  = 1:n;
stiff = P*pi^2/L*diag(tmp.^2)+E*I*pi^4/L^3*diag(tmp.^4);
% mu = (b,h,L,rho,E,P)
DM = {rho*h*L*eye(n)
    rho*b*L*eye(n)
    rho*A*eye(n)
    A*L*eye(n)
    zeros(n)
    zeros(n)};
DK = {E*dIdb*pi^4/L^3*diag(tmp.^4)
    E*dIdh*pi^4/L^3*diag(tmp.^4)
    -P*pi^2/L^2*diag(tmp.^2)-3*E*I*pi^4/L^4*diag(tmp.^4)
    zeros(n)
    I*pi^4/L^3*diag(tmp.^4)
    pi^2/L*diag(tmp.^2)};

tmpfun = @(input) zeros(n,1);
fnl{1} = tmpfun;
fnl{2} = @(input) cubic_nonlinearity(input,E,b,h,L);


dfnl2  = {tmpfun,tmpfun,tmpfun,tmpfun,tmpfun,tmpfun};

df3db = @(input) dcubic_db(input,E,b,h,L);
df3dh = @(input) dcubic_dh(input,E,b,h,L);
df3dL = @(input) dcubic_dL(input,E,b,h,L);
df3drho  = tmpfun;
df3dE = @(input) dcubic_dE(input,E,b,h,L);
df3dP = tmpfun;
dfnl3 = {df3db,df3dh,df3dL,df3drho,df3dE,df3dP};

end

function y = cubic_nonlinearity(x,E,b,h,L)
A  = b*h;
v1 = x{1};
v2 = x{2};
v3 = x{3};
n  = numel(v1);
y  = zeros(n,1);
seq = (1:n)';
tmp = sum(v2.*v3.*seq.^2);
for i=1:n
    y(i) = tmp*v1(i)*i^2;
end
y = E*A*pi^4/(4*L^3)*y;
end

function y = dcubic_db(x,E,b,h,L)
dAdb = h;
v1 = x{1};
v2 = x{2};
v3 = x{3};
n  = numel(v1);
y  = zeros(n,1);
seq = (1:n)';
tmp = sum(v2.*v3.*seq.^2);
for i=1:n
    y(i) = tmp*v1(i)*i^2;
end
y = E*dAdb*pi^4/(4*L^3)*y;
end

function y = dcubic_dh(x,E,b,h,L)
dAdh = b;
v1 = x{1};
v2 = x{2};
v3 = x{3};
n  = numel(v1);
y  = zeros(n,1);
seq = (1:n)';
tmp = sum(v2.*v3.*seq.^2);
for i=1:n
    y(i) = tmp*v1(i)*i^2;
end
y = E*dAdh*pi^4/(4*L^3)*y;
end

function y = dcubic_dL(x,E,b,h,L)
A  = b*h;
v1 = x{1};
v2 = x{2};
v3 = x{3};
n  = numel(v1);
y  = zeros(n,1);
seq = (1:n)';
tmp = sum(v2.*v3.*seq.^2);
for i=1:n
    y(i) = tmp*v1(i)*i^2;
end
y = -3*E*A*pi^4/(4*L^4)*y;
end

function y = dcubic_dE(x,E,b,h,L)
A  = b*h;
v1 = x{1};
v2 = x{2};
v3 = x{3};
n  = numel(v1);
y  = zeros(n,1);
seq = (1:n)';
tmp = sum(v2.*v3.*seq.^2);
for i=1:n
    y(i) = tmp*v1(i)*i^2;
end
y = A*pi^4/(4*L^3)*y;
end

