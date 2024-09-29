function [W,R] = explicit_whisker_InResonance(obj)
% Invariant manifold in the autonomous system limit

% setup
W  = struct();
R  = struct();
Ap = obj.System.A;
Bp = obj.System.B;
Q  = obj.System.F{1};
C  = obj.System.F{2};
midx = obj.System.N;
if ~obj.System.nl_damp && obj.System.order==2 
    midx = midx/2;  % 2nd order fct handle only takes displacement as input
end
ispropDamp = obj.System.Options.RayleighDamping; % damping type
Lambda_E = obj.E.spectrum;
assert(numel(Lambda_E)==4,'current implementation only supports four-dimensional SSMs');
idx = 1;
if imag(Lambda_E(1))<0; idx=2; end
Lambda_E = Lambda_E(idx);
vE  = obj.E.basis(:,idx); vEbar = conj(vE);
uE  = obj.E.adjointBasis(:,idx);

switch obj.System.order
    case 1
        %% first order computation
        % Linear terms
        W10 = vE;
        
        % quadratic order terms
        solver = obj.Options.solver;
        W20 = solveinveq(2*Lambda_E*Bp-Ap,Q({vE(1:midx);vE(1:midx)}),solver);
        W11 = solveinveq(2*real(Lambda_E)*Bp-Ap,...
            Q({vE(1:midx);vEbar(1:midx)})+Q({vEbar(1:midx);vE(1:midx)}),solver);
        
        % cubic order terms
        W30 = solveinveq(3*Lambda_E*Bp-Ap,Q({vE(1:midx);W20(1:midx)})+...
            Q({W20(1:midx);vE(1:midx)})+C({vE(1:midx);vE(1:midx);vE(1:midx)}),solver);
        tmp = Q({vEbar(1:midx);W20(1:midx)})+Q({W20(1:midx);vEbar(1:midx)})+...
            Q({vE(1:midx);W11(1:midx)})+Q({W11(1:midx);vE(1:midx)})+...
            C({vEbar(1:midx);vE(1:midx);vE(1:midx)})+...
            C({vE(1:midx);vEbar(1:midx);vE(1:midx)})+...
            C({vE(1:midx);vE(1:midx);vEbar(1:midx)});
        r21 = uE'*tmp;
        W21 = solveinveq((2*Lambda_E+conj(Lambda_E))*Bp-Ap,tmp-Bp*vE*r21,solver);

  
end

W.W10 = W10;
W.W20 = W20;
W.W11 = W11;
W.W30 = W30;
W.W21 = W21;
R.r21 = r21;
R.Lambda_E = Lambda_E;
end

