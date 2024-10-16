function [W,R] = explicit_whisker_InResonance(obj)
% Invariant manifold in the autonomous system limit

% setup
W  = struct();
R  = struct();
Ap = obj.System.A;
Bp = obj.System.B;
Q  = obj.System.F_semi{2};
C  = obj.System.F_semi{3};
midx = obj.System.N;
if ~obj.System.nl_damp && obj.System.order==2 
    midx = midx/2;  % 2nd order fct handle only takes displacement as input
end
ispropDamp = obj.System.Options.RayleighDamping; % damping type
Lambda_E = obj.E.spectrum;
assert(numel(Lambda_E)==4,'current implementation only supports four-dimensional SSMs');
idx1 = 1; idx2 = 3;
if imag(Lambda_E(idx1))<0; idx1=2; end
if imag(Lambda_E(idx2))<0; idx2=4; end
% if imag(Lambda_E(idx1))>imag(Lambda_E(idx2))
%     temp = idx1;
%     idx1 = idx2;
%     idx2 = temp;
% end
Lambda1_E = Lambda_E(idx1); Lambda1bar_E = conj(Lambda1_E);
Lambda2_E = Lambda_E(idx2); Lambda2bar_E = conj(Lambda2_E);
vE1  = obj.E.basis(:,idx1); vE1bar = conj(vE1);
vE2  = obj.E.basis(:,idx2); vE2bar = conj(vE2);
uE1  = obj.E.adjointBasis(:,idx1); uE1bar = conj(uE1);
uE2  = obj.E.adjointBasis(:,idx2); 

switch obj.System.order
    case 1
        %% first order computation
        % Linear terms
        W1000 = vE1;
        W0010 = vE2;
        
        % quadratic order terms
        solver = obj.Options.solver;
        r2000 = uE2'*Q({vE1(1:midx);vE1(1:midx)});
        gamma1 = r2000;
        W2000 = solveinveq(2*Lambda1_E*Bp-Ap,Q({vE1(1:midx);vE1(1:midx)})-gamma1*Bp*vE2,solver);
%         W0200 = conj(W2000);

        W0020 = solveinveq(2*Lambda2_E*Bp-Ap,Q({vE2(1:midx);vE2(1:midx)}),solver);
%         W0002 = conj(W0020);

        W1100 = solveinveq(2*real(Lambda1_E)*Bp-Ap,...
            Q({vE1(1:midx);vE1bar(1:midx)})+Q({vE1bar(1:midx);vE1(1:midx)}),solver);

        W1010 = solveinveq((Lambda1_E+Lambda2_E)*Bp-Ap,...
            Q({vE1(1:midx);vE2(1:midx)})+Q({vE2(1:midx);vE1(1:midx)}),solver);
%         W0101 = conj(W1010);

        r1001 = uE1bar'*(Q({vE1(1:midx);vE2bar(1:midx)})+Q({vE2bar(1:midx);vE1(1:midx)}));
        gamma2 = r1001;
        W1001 = solveinveq((Lambda1_E+Lambda2bar_E)*Bp-Ap,...
            Q({vE1(1:midx);vE2bar(1:midx)})+Q({vE2bar(1:midx);vE1(1:midx)})-gamma2*Bp*vE1bar,solver);
        W0110 = conj(W1001);

        W0011 = solveinveq(2*real(Lambda2_E)*Bp-Ap,...
            Q({vE2(1:midx);vE2bar(1:midx)})+Q({vE2bar(1:midx);vE2(1:midx)}),solver);


        % cubic order terms
        W3000 = solveinveq(3*Lambda1_E*Bp-Ap,Q({vE1(1:midx);W2000(1:midx)})+...
            Q({W2000(1:midx);vE1(1:midx)})+C({vE1(1:midx);vE1(1:midx);vE1(1:midx)})-gamma1*Bp*W1010,solver);

        W0030 = solveinveq(3*Lambda2_E*Bp-Ap,Q({vE2(1:midx);W0020(1:midx)})+...
            Q({W0020(1:midx);vE2(1:midx)})+C({vE2(1:midx);vE2(1:midx);vE2(1:midx)}),solver);

        tmp1 = Q({vE1(1:midx);W1100(1:midx)})+Q({W1100(1:midx);vE1(1:midx)})+...
            Q({vE1bar(1:midx);W2000(1:midx)})+Q({W2000(1:midx);vE1bar(1:midx)})+...
            C({vE1(1:midx);vE1(1:midx);vE1bar(1:midx)})+...
            C({vE1(1:midx);vE1bar(1:midx);vE1(1:midx)})+...
            C({vE1bar(1:midx);vE1(1:midx);vE1(1:midx)})-gamma1*Bp*W0110;
        r2100 = uE1'*tmp1;
        gamma3 = r2100;
        W2100 = solveinveq((2*Lambda1_E+Lambda1bar_E)*Bp-Ap,tmp1-gamma3*Bp*vE1,solver);

        W2010 = solveinveq((2*Lambda1_E+Lambda2_E)*Bp-Ap,Q({vE1(1:midx);W1010(1:midx)})+...
            Q({W1010(1:midx);vE1(1:midx)})+Q({W2000(1:midx);vE2(1:midx)})+...
            Q({vE2(1:midx);W2000(1:midx)})+C({vE1(1:midx);vE1(1:midx);vE2(1:midx)})+...
            C({vE1(1:midx);vE2(1:midx);vE1(1:midx)})+...
            C({vE2(1:midx);vE1(1:midx);vE1(1:midx)})-2*gamma1*Bp*W0020,solver);

        W2001 = solveinveq((2*Lambda1_E+Lambda2bar_E)*Bp-Ap,Q({vE1(1:midx);W1001(1:midx)})+...
            Q({W1001(1:midx);vE1(1:midx)})+Q({W2000(1:midx);vE2bar(1:midx)})+...
            Q({vE2bar(1:midx);W2000(1:midx)})+C({vE1(1:midx);vE1(1:midx);vE2bar(1:midx)})+...
            C({vE1(1:midx);vE2bar(1:midx);vE1(1:midx)})+...
            C({vE2bar(1:midx);vE1(1:midx);vE1(1:midx)})-gamma1*Bp*W0011-gamma2*Bp*W1100,solver);

        W1020 = solveinveq((Lambda1_E+2*Lambda2_E)*Bp-Ap,Q({vE2(1:midx);W1010(1:midx)})+...
            Q({W1010(1:midx);vE2(1:midx)})+Q({W0020(1:midx);vE1(1:midx)})+...
            Q({vE1(1:midx);W0020(1:midx)})+C({vE2(1:midx);vE2(1:midx);vE1(1:midx)})+...
            C({vE2(1:midx);vE1(1:midx);vE2(1:midx)})+...
            C({vE1(1:midx);vE2(1:midx);vE2(1:midx)}),solver);

        W0120 = solveinveq((Lambda1bar_E+2*Lambda2_E)*Bp-Ap,Q({vE2(1:midx);W0110(1:midx)})+...
            Q({W0110(1:midx);vE2(1:midx)})+Q({W0020(1:midx);vE1bar(1:midx)})+...
            Q({vE1bar(1:midx);W0020(1:midx)})+C({vE2(1:midx);vE2(1:midx);vE1bar(1:midx)})+...
            C({vE2(1:midx);vE1bar(1:midx);vE2(1:midx)})+...
            C({vE1bar(1:midx);vE2(1:midx);vE2(1:midx)})-conj(gamma2)*Bp*W1010,solver);

        tmp2  = Q({vE2(1:midx);W0011(1:midx)})+...
            Q({W0011(1:midx);vE2(1:midx)})+Q({W0020(1:midx);vE2bar(1:midx)})+...
            Q({vE2bar(1:midx);W0020(1:midx)})+C({vE2(1:midx);vE2(1:midx);vE2bar(1:midx)})+...
            C({vE2(1:midx);vE2bar(1:midx);vE2(1:midx)})+...
            C({vE2bar(1:midx);vE2(1:midx);vE2(1:midx)});
        r0021 = uE2'*tmp2;
        gamma4 = r0021;
        W0021 = solveinveq((2*Lambda2_E+Lambda2bar_E)*Bp-Ap,tmp2-gamma4*Bp*vE2,solver);

        tmp3 = Q({vE1(1:midx);W0110(1:midx)})+Q({vE1bar(1:midx);W1010(1:midx)})...
            +Q({vE2(1:midx);W1100(1:midx)})+Q({W0110(1:midx);vE1(1:midx)})...
            +Q({W1010(1:midx);vE1bar(1:midx)})+Q({W1100(1:midx);vE2(1:midx)})...
            +C({vE1(1:midx);vE1bar(1:midx);vE2(1:midx)})...
            +C({vE1(1:midx);vE2(1:midx);vE1bar(1:midx)})...
            +C({vE1bar(1:midx);vE1(1:midx);vE2(1:midx)})...
            +C({vE1bar(1:midx);vE2(1:midx);vE1(1:midx)})...
            +C({vE2(1:midx);vE1(1:midx);vE1bar(1:midx)})...
            +C({vE2(1:midx);vE1bar(1:midx);vE1(1:midx)})-2*conj(gamma2)*Bp*W2000;
        r1110 = uE2'*tmp3;
        gamma5 = r1110;
        W1110 = solveinveq((Lambda1_E+Lambda1bar_E+Lambda2_E)*Bp-Ap,tmp3-gamma5*Bp*vE2,solver);

        tmp4 = Q({vE1(1:midx);W0011(1:midx)})+Q({vE2(1:midx);W1001(1:midx)})...
            +Q({vE2bar(1:midx);W1010(1:midx)})+Q({W0011(1:midx);vE1(1:midx)})...
            +Q({W1001(1:midx);vE2(1:midx)})+Q({W1010(1:midx);vE2bar(1:midx)})...
            +C({vE1(1:midx);vE2(1:midx);vE2bar(1:midx)})...
            +C({vE1(1:midx);vE2bar(1:midx);vE2(1:midx)})...
            +C({vE2(1:midx);vE1(1:midx);vE2bar(1:midx)})...
            +C({vE2(1:midx);vE2bar(1:midx);vE1(1:midx)})...
            +C({vE2bar(1:midx);vE1(1:midx);vE2(1:midx)})...
            +C({vE2bar(1:midx);vE2(1:midx);vE1(1:midx)})-gamma2*Bp*W0110;
        r1011 = uE1'*tmp4;
        gamma6 = r1011;
        W1011 = solveinveq((Lambda1_E+Lambda2_E+Lambda2bar_E)*Bp-Ap,tmp4-gamma6*Bp*vE1,solver);






end
% Linear terms
W.W1000 = W1000;
W.W0010 = W0010;

% Quadratic order terms
W.W2000 = W2000;
W.W0020 = W0020;
W.W1100 = W1100;
W.W1010 = W1010;
W.W1001 = W1001;
W.W0011 = W0011;
R.r2000 = r2000;
R.r1001 = r1001;

% Cubic order terms
W.W3000 = W3000;
W.W0030 = W0030;
W.W2100 = W2100;
W.W2010 = W2010;
W.W2001 = W2001;
W.W1020 = W1020;
W.W0120 = W0120;
W.W0021 = W0021;
W.W1110 = W1110;
W.W1011 = W1011;
R.r2100 = r2100;
R.r0021 = r0021;
R.r1110 = r1110;
R.r1011 = r1011;

R.Lambda_E = Lambda_E;
end

