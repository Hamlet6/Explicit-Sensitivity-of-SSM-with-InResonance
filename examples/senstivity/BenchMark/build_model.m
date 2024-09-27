function [mass,damp,gyro,stiff,fnl] = build_model(a,b,c,d)

n = 2;
mass = eye(2);
damp = [0.01 0;0 0.03];
gyro = a*[0 1;-1 0];
stiff = [1 0;0 2];

subs2 = [1 1 2
    1 3 3
    2 1 1];
vals2 = [b; c; b];
subs3 = [1 1 1 1
    1 1 2 2
    2 2 2 2
    2 3 3 4];
vals3 = [1; 1; 1.6; d];

f2 = sptensor(subs2, vals2, [n,2*n,2*n]); 
f3 = sptensor(subs3, vals3, [n,2*n,2*n,2*n]); 
fnl = {f2,f3};

end