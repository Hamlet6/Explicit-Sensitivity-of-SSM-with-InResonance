function [mass,damp,stiff,fnl,fext] = build_model(mu,f,n)

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
cubic_coeff = zeros(n,n,n,n);
quad_coeff1  = zeros(n,n,n);
quad_coeff2  = zeros(n,n,n);
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

for i = 1: n
    for j = 1:n
        for k = 1:n
            for l= 1:n
                cubic_coeff(i,j,k,l) = 0.5*B(i,j)*B(k,l);
            end
        end
    end
end

for i = 1: n
    for j = 1:n
        for m = 1:n
            quad_coeff1(i,j,m) = B(i,j)*d(m);
        end
    end
end

for i = 1: n
    for k = 1:n
        for l = 1:n
            quad_coeff2(i,k,l) = d(i)*B(k,l);
        end
    end
end

tmp2  = sptensor(quad_coeff1);
tmp3  = sptensor(quad_coeff2);
subs2 = tmp2.subs;
subs3 = tmp3.subs;
subs23 = [subs2; subs3];
vals23 = [tmp2.vals; 0.5*tmp3.vals];
fnl23 = sptensor(subs23, vals23, [n,n,n]);

fnl1 = sptensor(cubic_coeff);

mass = M;
damp = C;
stiff = K+H;
fnl = {fnl23,fnl1};
fext = F;

end