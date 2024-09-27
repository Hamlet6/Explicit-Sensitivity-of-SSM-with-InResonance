function [mass,damp,gyro,stiff,fnl,DM,DC,DK,dfnl2,dfnl3] = build_model_semi(n,alpha,gamma,k1,kf)

% A = 0.04*0.03; % m^2
% I = 0.04*0.03^3/12;
% rho = 7680;
% E = 30e9;
% % eta = 1e-5*E;
% eta = 1e-4*E;
% L = 1;
% P = 67.5e3;
% 
% kf = sqrt(E*I/(P*L^2));
% k1 = sqrt(E*A/P);
% alpha = I*eta/(L^3*sqrt(rho*A*P));
mass  = eye(n);
damp  = zeros(n);
gyro  = zeros(n);
stiff = zeros(n);
fext  = zeros(n,1);
cubic_coeff = zeros(n,n,n,n);
for i=1:n
    damp(i,i)  = alpha*(i*pi)^4;
    stiff(i,i) = kf^2*(i*pi)^4-(gamma^2-1)*(i*pi)^2;
    fext(i)    = (1-(-1)^i)/(i*pi);
    for j=1:n
        if j~=i
            gyro(i,j) = 4*gamma*i*j/(i^2-j^2)*(1-(-1)^(i+j));
        end
        cubic_coeff(i,j,j,i) = i^2*j^2;
    end
end

% mu = (alpha,gamma,k1,kf)
DM = {zeros(n)
    zeros(n)
    zeros(n)
    zeros(n)};
DC = {damp/alpha
    gyro/gamma
    zeros(n)
    zeros(n)};
tmp = (1:n)*pi;
DK = {zeros(n)
    -2*gamma*diag(tmp.^2)
    zeros(n)
    2*kf*diag(tmp.^4)};

disp('the first four eigenvalues for undamped system');
lamd = eigs([zeros(n) eye(n);-mass\stiff -mass\(gyro)],4,'smallestabs')

tmpfun = @(input) zeros(n,1);
fnl{1} = tmpfun;
fnl{2} = @(input) cubic_nonlinearity(input,k1,alpha,kf);

dfnl2  = {tmpfun,tmpfun,tmpfun,tmpfun};

df3dalpha = @(input) dcubic_dalpha(input,k1,alpha,kf);
df3dgamma = tmpfun;
df3dk1  = @(input) dcubic_dk1(input,k1,alpha,kf);
df3dkf  = @(input) dcubic_dkf(input,k1,alpha,kf);
dfnl3 = {df3dalpha,df3dgamma,df3dk1,df3dkf};

end

function y = cubic_nonlinearity(z,k1,alpha,kf)
z1 = z{1}; 
z2 = z{2};
z3 = z{3};
n  = numel(z1)/2;
x1 = z1(1:n); 
x2 = z2(1:n); v2 = z2(n+1:2*n);
x3 = z3(1:n); 
coef1 = 0.25*k1^2*pi^4;
coef2 = 0.5*alpha*k1^2/kf^2*pi^4;
y   = zeros(n,1);
seq = (1:n)';
tmp1 = sum(x1.*x2.*seq.^2);
tmp2 = sum(x1.*v2.*seq.^2);
for i=1:n
    y(i) = coef1*tmp1*x3(i)*i^2+coef2*tmp2*x3(i)*i^2;
end
end

function y = dcubic_dalpha(z,k1,alpha,kf)
z1 = z{1}; 
z2 = z{2};
z3 = z{3};
n  = numel(z1)/2;
x1 = z1(1:n); 
v2 = z2(n+1:2*n);
x3 = z3(1:n); 
coef2 = 0.5*k1^2/kf^2*pi^4;
y   = zeros(n,1);
seq = (1:n)';
tmp2 = sum(x1.*v2.*seq.^2);
for i=1:n
    y(i) = coef2*tmp2*x3(i)*i^2;
end
end

function y = dcubic_dk1(z,k1,alpha,kf)
z1 = z{1}; 
z2 = z{2};
z3 = z{3};
n  = numel(z1)/2;
x1 = z1(1:n); 
x2 = z2(1:n); v2 = z2(n+1:2*n);
x3 = z3(1:n); 
coef1 = 0.5*k1*pi^4;
coef2 = alpha*k1/kf^2*pi^4;
y   = zeros(n,1);
seq = (1:n)';
tmp1 = sum(x1.*x2.*seq.^2);
tmp2 = sum(x1.*v2.*seq.^2);
for i=1:n
    y(i) = coef1*tmp1*x3(i)*i^2+coef2*tmp2*x3(i)*i^2;
end
end

function y = dcubic_dkf(z,k1,alpha,kf)
z1 = z{1}; 
z2 = z{2};
z3 = z{3};
n  = numel(z1)/2;
x1 = z1(1:n); 
v2 = z2(n+1:2*n);
x3 = z3(1:n); 
coef2 = -alpha*k1^2/kf^3*pi^4;
y   = zeros(n,1);
seq = (1:n)';
tmp2 = sum(x1.*v2.*seq.^2);
for i=1:n
    y(i) = coef2*tmp2*x3(i)*i^2;
end
end