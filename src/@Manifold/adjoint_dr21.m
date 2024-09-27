function Dr21 = adjoint_dr21(obj,W20,W11,W30,W21,r21)
% Computation of senstivity of the cubic term coefficient using adjoint method. 
% The current implementation is restricted to second order system with symmetric
% mass and stiffness matrices along with propotional damping. 

%% setup 
% check system setup
is2ndorder = obj.System.order==2; % second order system
ispropDamp = obj.System.Options.RayleighDamping; % damping type
% spectrum
Lambda_E = obj.E.spectrum;
assert(numel(Lambda_E)==2,'current implementation only supports two-dimensional SSMs');
idx = 1;
if imag(Lambda_E(1))<0; idx=2; end
Lambda_E = Lambda_E(idx); midx = obj.System.n;
vE  = obj.E.basis(:,idx); uE = obj.E.adjointBasis(:,idx);
phi = vE(1:midx);
psi = conj(uE(1:midx));
% second-order system
M   = obj.System.M;
K   = obj.System.K;
C   = obj.System.C;
f2  = obj.System.fnl_semi{1};
f2  = @(x1,x2) f2({x1;x2});
f3  = obj.System.fnl_semi{2};
f3  = @(x1,x2,x3) f3({x1;x2;x3});

% derivative w.r.t parameters
DM = obj.System.DM;
DK = obj.System.DK;
DC = obj.System.DC;
Dfnl2 = obj.System.Dfnl2_semi;
Dfnl3 = obj.System.Dfnl3_semi;
nps = numel(DM); % number of parameters
Df2 = cell(nps,1);
Df3 = cell(nps,1);
for k = 1:nps
    Df2{k} = @(x1,x2) Dfnl2{k}({x1;x2});
    Df3{k} = @(x1,x2,x3) Dfnl3{k}({x1;x2;x3});
end

%% direct computation of adjoints
lamd21 = psi;       % \delta f21

%% inverse of n-by-n matrix
RLambda_E = real(Lambda_E);
A20 = 4*Lambda_E^2*M+2*Lambda_E*C+K;
A11 = 4*(RLambda_E)^2*M+2*RLambda_E*C+K;
A11hat = blkdiag(A11,A11);
[F22_vE,F23_vE] = F2func(vE,f2);
b11 = (F23_vE+F22_vE).'*lamd21;
solver = obj.Options.solver;
% lamd20 = solveinveq(A20.',bhat,solver);    % \delta W201
lamd11 = solveinveq(A11hat.',b11,solver);    % \delta W111
In = eye(midx);
F11_phibar = 0.5*[In; 2*RLambda_E*In]*(F22_vE+F23_vE)*[In; conj(Lambda_E)*In];

[F22_W20,F23_W20] = F2func(W20,f2);
[F32_vE,F33_vE,F34_vE] = F3func(vE,f3);
F3_vbar = F32_vE+F33_vE+F34_vE;
f21_phibar = (F22_W20+F23_W20+F3_vbar)*[In; conj(Lambda_E)*In];
tmp = -f21_phibar.'*lamd21+2*F11_phibar.'*lamd11;
%% inverse of (n+1)-by-(n+1) matrix

%% calculate gradient
Dr21 = zeros(nps,1); 
for k=1:nps
    temp1 = -lamd21.'*(Df2{k}(phi,W201)+Df2{k}(W201,phi)+Df2{k}(phi,W111)+Df2{k}(W111,phi)+3*Df3{k}(phi,phi,phi))+...
        (lamd20+2*lamd11).'*Df2{k}(phi,phi);
    temp2 = lamd20.'*((4*Lambda_E^2+2*Lambda_E*al)*DM{k}+(2*Lambda_E*be+1)*DK{k})*W201;
    temp3 = lamd11.'*((4*RLambda_E^2+2*RLambda_E*al)*DM{k}+(2*RLambda_E*be+1)*DK{k})*W111;
    temp4 = lamd_phi.'*(DK{k}-om^2*DM{k})*phi+lamd_norm*phi'*DM{k}*phi;
    Dr21(k) = temp1+temp2+temp3+temp4;
end

end


function [F22,F23] = F2func(x,f2)
% F2FUNC This function returns two matrices
n = numel(x);
E = eye(n);
% the following implementation is not optimized. One can consider sparse
% matrix along with parallel computation.
F22 = zeros(n/2,n);
F23 = zeros(n/2,n);
for k=1:n
    ek = E(:,k);
    F22(:,k) = f2(ek,x);
    F23(:,k) = f2(x,ek);
end
end



function [F32,F33,F34] = F3func(x,f3)
% F3FUNC This function returns two matrices
n = numel(x);
E = eye(n);
% the following implementation is not optimized. One can consider sparse
% matrix along with parallel computation.
F32 = zeros(n/2,n); F33 = zeros(n/2,n); F34 = zeros(n/2,n);
for k=1:n
    ek = E(:,k);
    F32(:,k) = f3(ek,x,x);
    F33(:,k) = f3(x,ek,x);
    F34(:,k) = f3(x,x,ek);
end
end