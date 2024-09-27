function [mass,damp,gyro,stiff,fnl,DM,DC,DK,dfnl2,dfnl3] = build_model_semi(a,b,c,d,isgyro)

n = 2;
mass = eye(2);
damp = [0.01 0;0 0.03];
stiff = [1 0;0 2];

if isgyro
    gyro = a*[0 1;-1 0];
    % mu = (a,b,c,d)
    DM = {zeros(n)
        zeros(n)
        zeros(n)
        zeros(n)};
    DC = {gyro/a
        zeros(n)
        zeros(n)
        zeros(n)};
    DK = DM;
    
    % nonlinear function handles and their derivatives
    fnl{1} = @(input) quadratic_nonlinearity(input,b,c);
    fnl{2} = @(input) cubic_nonlinearity(input,d);
    
    tmpfun  = @(input) zeros(n,1);
    dfnl2db = @(input) dquadratic_db(input,b,c);
    dfnl2dc = @(input) dquadratic_dc(input,b,c);
    dfnl2   = {tmpfun,dfnl2db,dfnl2dc,tmpfun};
    dfnl3dd = @(input) dcubic_dd(input,d);
    dfnl3   = {tmpfun,tmpfun,tmpfun,dfnl3dd};
else
    % no gyroscopic forces
    gyro = zeros(2);
    % mu = (b,c,d)
    DM = {zeros(n)
        zeros(n)
        zeros(n)};
    DC = {zeros(n)
        zeros(n)
        zeros(n)};
    DK = DM;
    
    % nonlinear function handles and their derivatives
    fnl{1} = @(input) quadratic_nonlinearity(input,b,c);
    fnl{2} = @(input) cubic_nonlinearity(input,d);
    
    tmpfun  = @(input) zeros(n,1);
    dfnl2db = @(input) dquadratic_db(input,b,c);
    dfnl2dc = @(input) dquadratic_dc(input,b,c);
    dfnl2   = {dfnl2db,dfnl2dc,tmpfun};
    dfnl3dd = @(input) dcubic_dd(input,d);
    dfnl3   = {tmpfun,tmpfun,dfnl3dd};
end

end

function y = quadratic_nonlinearity(z,b,c)
z1 = z{1}; 
z2 = z{2};
y = [b*z1(1)*z2(2)+c*z1(3)*z2(3)
    b*z1(1)*z2(1)];
end

function y = dquadratic_db(z,b,c)
z1 = z{1}; 
z2 = z{2};
y = [z1(1)*z2(2)
    z1(1)*z2(1)];
end

function y = dquadratic_dc(z,b,c)
z1 = z{1}; 
z2 = z{2};
y = [z1(3)*z2(3)
    0];
end

function y = cubic_nonlinearity(z,d)
z1 = z{1}; 
z2 = z{2};
z3 = z{3};
y = [z1(1)*z2(1)*z3(1)+z1(1)*z2(2)*z3(2)
    1.6*z1(2)*z2(2)*z3(2)+d*z1(3)*z2(3)*z3(4)];
end

function y = dcubic_dd(z,d)
z1 = z{1}; 
z2 = z{2};
z3 = z{3};
y = [0
    z1(3)*z2(3)*z3(4)];
end