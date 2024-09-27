function [mass,damp,stiff,fnl,fext] = build_model_semi(c1,c2,b1,b2,f1,f2)

n = 2;
mass = eye(n,n);
damp = [c1, 0;
     0, c2];
stiff = [1 0;0 2];
subs2 = [1 1 2
    2 1 1];
vals2 = [b1 b2]';
F2 = sptensor(subs2, vals2, [n,n,n]);
fnl = {F2};

fnl{1} = @(input) quadratic_nonlinearity(input,b1,b2);
fnl{2} = @(input) cubic_nonlinearity(input);

fext = [f1;f2];

end

function y = quadratic_nonlinearity(x,b1,b2)
v1 = x{1};
v2 = x{2};
y = [b1*v1(1)*v2(2); b2*v1(1)*v2(1)];
end

function y = cubic_nonlinearity(x)
v1 = x{1};
v2 = x{2};
v3 = x{3};
y = [v1(1)*v2(1)*v3(1)+v1(1)*v2(2)*v3(2)
    0.3*v1(1)*v2(2)*v3(2)+1.6*v1(2)*v2(2)*v3(2)];
end