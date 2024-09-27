function Ln = internalforce(model, u, Null, Nullf, ud)
u0 = Null*u + ud;
% x=u;
model.sol('sol1').setU(u0);
model.sol('sol1').setPVals(0);
model.sol('sol1').createSolution;
ML = mphmatrix(model, 'sol1', 'Out', {'L','K'}, 'initmethod','sol', 'initsol', 'sol1');
% Kn = ML.K;
Ln = Nullf'*(ML.L-ML.K*ud);


end