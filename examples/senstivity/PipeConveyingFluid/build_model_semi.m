function [mass,damp,stiff,fnl,fext,DM,DC,DK,dfnl2,dfnl3,dfext] = build_model_semi(n,u,beta,alpha,BC)

l = 1;
disp('Getting linear and nonlinearity coefficients')
fileName = [BC,'_tensors_',num2str(n),'_',num2str(1),'.mat'];
try 
    load(fileName,'a_ijkl','b_ijkl','c_ijkl','g_ijkl','a','b','delta','delta1','intphi');   
    disp('Loaded coefficients from storage');
catch
    disp('Calculating coefficients');
    [a,b,delta,delta1,a_ijkl,b_ijkl,c_ijkl,g_ijkl,intphi] = cal_parameters(n,l,BC);
    disp('Saving coefficients');
    save(fileName,'a_ijkl','b_ijkl','c_ijkl','g_ijkl','a','b','delta','delta1','intphi','-v7.3')
end

mass  = delta;
damp  = alpha*delta1+2*u*sqrt(beta)*a;
stiff = delta1+u^2*b;
fext  = intphi;

% mu = (u,beta,alpha)
DM = {zeros(n)
    zeros(n)
    zeros(n)};
DC = {2*sqrt(beta)*a
    u/sqrt(beta)*a
    delta1};
DK = {2*u*b
    zeros(n)
    zeros(n)};
dfext = {zeros(n,1),zeros(n,1),zeros(n,1)};

% nonlinear terms
alpha_ijkl = u^2*c_ijkl+a_ijkl;
beta_ijkl  = 2*u*sqrt(beta)*b_ijkl;
gamma_ijkl = g_ijkl;
tmpfun = @(input) zeros(n,1);
fnl{1} = tmpfun;
fnl{2} = @(input) cubic_nonlinearity(input,alpha_ijkl,beta_ijkl,gamma_ijkl);

dfnl2  = {tmpfun,tmpfun,tmpfun};

df3du      = @(input) dcubic_du(input,2*u*c_ijkl,2*sqrt(beta)*b_ijkl);
df3dbeta   = @(input) dcubic_dbeta(input,u/sqrt(beta)*b_ijkl);
df3dalpha  = tmpfun;
dfnl3      = {df3du,df3dbeta,df3dalpha};
        
end

function y = cubic_nonlinearity(z,alpha_ijkl,beta_ijkl,gamma_ijkl)
z1 = z{1}; 
z2 = z{2};
z3 = z{3};
n  = numel(z1)/2;
x1 = z1(1:n); 
x2 = z2(1:n); v2 = z2(n+1:2*n);
x3 = z3(1:n); v3 = z3(n+1:2*n);
y   = zeros(n,1);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                y(i) = y(i)+alpha_ijkl(i,j,k,l)*x1(j)*x2(k)*x3(l)+...
                    beta_ijkl(i,j,k,l)*x1(j)*x2(k)*v3(l)+...
                    gamma_ijkl(i,j,k,l)*x1(j)*v2(k)*v3(l);
            end
        end
    end
end
end

function y = dcubic_du(z,alpha_ijkl,beta_ijkl)
z1 = z{1}; 
z2 = z{2};
z3 = z{3};
n  = numel(z1)/2;
x1 = z1(1:n); 
x2 = z2(1:n); 
x3 = z3(1:n); v3 = z3(n+1:2*n);
y   = zeros(n,1);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                y(i) = y(i)+alpha_ijkl(i,j,k,l)*x1(j)*x2(k)*x3(l)+...
                    beta_ijkl(i,j,k,l)*x1(j)*x2(k)*v3(l);
            end
        end
    end
end
end

function y = dcubic_dbeta(z,beta_ijkl)
z1 = z{1}; 
z2 = z{2};
z3 = z{3};
n  = numel(z1)/2;
x1 = z1(1:n); 
x2 = z2(1:n); 
v3 = z3(n+1:2*n);
y   = zeros(n,1);
for i=1:n
    for j=1:n
        for k=1:n
            for l=1:n
                y(i) = y(i)+beta_ijkl(i,j,k,l)*x1(j)*x2(k)*v3(l);
            end
        end
    end
end
end