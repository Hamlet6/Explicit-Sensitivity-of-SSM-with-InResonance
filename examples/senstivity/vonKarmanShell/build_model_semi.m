function [M,K,fnl,fext,DM,DK,dfnl2,dfnl3,dfext,outdof] = build_model_semi(nDiscretization,E,rho)

% build unit model
[M0,K0,fnl0,dfnl0,fext,outdof] = build_unit_model_semi(nDiscretization);

% recale matrices and tensors
M = M0*rho;
K = K0*E;
fnl = cell(1,2);
fnl{1}  = @(input) E*fnl0{1}(input);
fnl{2}  = @(input) E*fnl0{2}(input);

% derivatives w.r.t. rho and E
n = size(M0,1);
DM = {M0, zeros(n)};
DK = {zeros(n), K0};

tmp = @(input) zeros(n,1);
dfnl2dE   = @(input) dfnl0{1}(input);
dfnl2     = {tmp, dfnl2dE};
dfnl3dE   = @(input) dfnl0{2}(input);
dfnl3     = {tmp, dfnl3dE};
dfext     = {zeros(n,1), zeros(n,1)};

end
