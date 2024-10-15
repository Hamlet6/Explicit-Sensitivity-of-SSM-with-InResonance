function [mass,damp,stiff,fnl,fext] = build_model_semi(mu,f,n)

syms x 
F0 = f*dirac(x-0.25);   % excitation at x=0.25
a0=14.2365;%2.3788/14.2365
Z0=4*a0*x*(1-x);


M=zeros(n,n);
K=zeros(n,n);
H=zeros(n,n);
C=zeros(n,n);
B=zeros(n,n);
F=zeros(n,1);
d=zeros(n,1);
for i = 1:n
    for j = 1:n
        beita_i = i*pi;
        beita_j = j*pi;
        phi_i = sqrt(2)*sin(beita_i*x);
        phi_j = sqrt(2)*sin(beita_j*x);
        M(i,j) = int(phi_i*phi_j,0,1);
        K(i,j) = int(phi_i*diff(phi_j,4),0,1);
        H(i,j) = int(phi_i*diff(Z0,2),0,1)*int(phi_j*diff(Z0,2),0,1);
        C(i,j) = int(2*mu*phi_i*phi_j,0,1);
        B(i,j) = int(phi_i*diff(phi_j,2),0,1);
    end
    F(i) = int(phi_i*F0,0,1);
    d(i) = int(phi_i*diff(Z0,2),0,1);
end


mass = M;
damp = C;
stiff = K+H;

fnl{1} = @(input) quadratic_nonlinearity(input,B,d,n);
fnl{2} = @(input) cubic_nonlinearity(input,B,n);

fext = F;


end

function y = quadratic_nonlinearity(q,B,d,n)
q1 = q{1};
q2 = q{2};
y = 0.5*q1(1:n).'*B*q2(1:n)*d + d'*q1(1:n)*B*q2(1:n);
end

function y = cubic_nonlinearity(q,B,n)
q1 = q{1};
q2 = q{2};
q3 = q{3};
y = 0.5*q1(1:n).'*B*q2(1:n)*B*q3(1:n);
end