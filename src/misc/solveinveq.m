function [W] = solveinveq(C,b,solver)
% solves C*W = R for W with given solver
switch solver
    case 'lsqminnorm'
        W = lsqminnorm(C,b);
    case 'linsolve'
        W = linsolve(C,b);
    case 'backslash'
        W = C \b;
    case 'pinv'
        W = pinv(C)*b;
    case 'inv'
        W = inv(C)*b;
end
end