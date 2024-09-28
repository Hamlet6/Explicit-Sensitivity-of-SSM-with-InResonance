function [mass,damp,stiff,fnl,fext,DM,DK,dfnl2,dfnl3,dfext] = build_model_sens_direct(m,b1,b2,zeta)

n = 2;
mass = m*eye(n,n);
stiff = [1 0;0 2];
damp  = 2*zeta*sqrt(m)*eye(2);
DM = {eye(n),zeros(n),zeros(n)};
DK = {zeros(n),zeros(n),zeros(n)};

fnl{1} = @(input) quadratic_nonlinearity(input,b1);
fnl{2} = @(input) cubic_nonlinearity(input,b2);

fext = [m;b1*b2];

df2dm  = @(input) zeros(n,1);
df2db1 = @(input) dquad_db1(input);
df2db2 = @(input) zeros(n,1);
dfnl2  = {df2dm,df2db1,df2db2};

df3dm  = @(input) zeros(n,1);
df3db1 = @(input) zeros(n,1);
df3db2 = @(input) dcubic_db2(input);
dfnl3  = {df3dm,df3db1,df3db2};

dfextdm  = [1;0];
dfextdb1 = [0; b2];
dfextdb2 = [0; b1];
dfext    = {dfextdm,dfextdb1,dfextdb2};
end

function y = quadratic_nonlinearity(x,b1)
v1 = x{1};
v2 = x{2};
y = [b1*v1(1)*v2(2); v1(1)*v2(1)];
end

function y = cubic_nonlinearity(x,b2)
v1 = x{1};
v2 = x{2};
v3 = x{3};
y = [v1(1)*v2(1)*v3(1)+v1(1)*v2(2)*v3(2)
    b2*v1(1)*v2(2)*v3(2)+1.6*v1(2)*v2(2)*v3(2)];
end

function y = dquad_db1(x)
v1 = x{1};
v2 = x{2};
y = [v1(1)*v2(2); 0];
end

function y = dcubic_db2(x)
v1 = x{1};
v2 = x{2};
v3 = x{3};
y = [0
    v1(1)*v2(2)*v3(2)];
end

