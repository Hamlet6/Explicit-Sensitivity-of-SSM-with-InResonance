function [mass,damp,stiff,fnl,fext] = build_model(m,zeta,b1,b2)

n = 2;
mass = m*eye(n,n);
damp = 2*zeta*sqrt(m)*eye(n,n);
stiff = [1 0;0 2];
subs2 = [1 1 2
    2 1 1];
vals2 = [b1 1]';
subs3 = [1 1 1 1
    1 1 2 2
    2 1 2 2
    2 2 2 2];
vals3 = [1 1 b2 1.6]';
F2 = sptensor(subs2, vals2, [n,n,n]);
F3 = sptensor(subs3, vals3, [n,n,n,n]);
fnl = {F2,F3};
fext = [m;b1*b2];

end