function Dr21 = adjoint_sentivity(obj,W20,W11,W30,W21,r21)
% Computation of senstivity of the cubic term coefficient using adjoint method. 
% The current implementation is restricted to second order system with symmetric
% mass and stiffness matrices along with propotional damping. 

%% setup 
% check system setup
is2ndorder = obj.System.order==2; % second order system
ispropDamp = ~obj.System.nl_damp && obj.System.Options.RayleighDamping; % damping type
assert(is2ndorder&&ispropDamp, "current implementation assumes 2nd system with proportional damping");
% spectrum
Lambda_E = obj.E.spectrum;
assert(numel(Lambda_E)==2,'current implementation only supports two-dimensional SSMs');
idx = 1;
if imag(Lambda_E(1))<0; idx=2; end
Lambda_E = Lambda_E(idx); midx = obj.System.n;
vE  = obj.E.basis(:,idx);
phi = vE(1:midx);
assert(isreal(phi),'eigenvector is not real');

% second-order system
M   = obj.System.M;
K   = obj.System.K;
tmp = phi'*M*phi;
phi = phi/sqrt(tmp);
om  = sqrt(phi'*K*phi); lamdhat = om^2;
al  = obj.System.al(om);
be  = obj.System.be(om);
C   = al*M+be*K;
zet = phi'*C*phi/(2*om);
fprintf('\n damping ratio is %d\n', zet);
assert(abs(Lambda_E-(-zet*om+1i*om*sqrt(1-zet^2)))<1e-6*om,'eigenvalue does not match');
dal = obj.System.dal(om);
dbe = obj.System.dbe(om);
f2  = obj.System.fnl_semi{1};
f2  = @(x1,x2) f2({x1;x2});
f3  = obj.System.fnl_semi{2};
f3  = @(x1,x2,x3) f3({x1;x2;x3});

% derivative w.r.t parameters
DM = obj.System.DM;
DK = obj.System.DK;
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
lamd_gamma = -1;                    % \delta\gamma
W201 = W20(1:midx);
W111 = W11(1:midx);
f21  = f2(phi,W201)+f2(W201,phi)+f2(phi,W111)+f2(W111,phi)+3*f3(phi,phi,phi);
lamd_omega = lamd_gamma*phi'*f21/om;           % \delta\kappa
kapa = 1i/(2*om*sqrt(1-zet^2));
lamd21 = lamd_gamma*kapa*phi;       % \delta f21

%% inverse of n-by-n matrix
RLambda_E = real(Lambda_E);
A20 = 4*Lambda_E^2*M+2*Lambda_E*C+K;
A11 = 4*(RLambda_E)^2*M+2*RLambda_E*C+K;
[F22_phi,F23_phi] = F2func(phi,f2);
bhat = (F23_phi+F22_phi).'*lamd21;
solver = obj.Options.solver;
lamd20 = solveinveq(A20.',bhat,solver);    % \delta W201
lamd11 = solveinveq(A11.',bhat,solver);    % \delta W111

%% inverse of (n+1)-by-(n+1) matrix
[F22_W201,F23_W201] = F2func(W201,f2);
[F22_W111,F23_W111] = F2func(W111,f2);
[F32_phi,F33_phi,F34_phi] = F3func(phi,f3);
f21phi = F22_W201+F23_W201+F22_W111+F23_W111+3*(F32_phi+F33_phi+F34_phi);
bphi  = lamd_gamma*kapa*f21+f21phi.'*lamd21-(F22_phi+F23_phi).'*(lamd20+2*lamd11);
tmp1  = lamd20.'*(8*Lambda_E*M+2*C)*W201*(-zet+1i*sqrt(1-zet^2))-zet*transpose(lamd11)*(8*RLambda_E*M+2*C)*W111;
tmp2  = dal*M+dbe*K;
tmp3  = 2*Lambda_E*transpose(lamd20)*tmp2*W201+2*RLambda_E*transpose(lamd11)*tmp2*W111+lamd_omega*kapa;
bnorm = tmp1+tmp3;
coefmat = [2*om*phi'*M 0; K-om^2*M 2*M*phi];
soltmp  = solveinveq(coefmat, [bnorm;bphi],solver);
lamd_phi  = soltmp(1:midx);
lamd_norm = soltmp(end);

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
F22 = zeros(n);
F23 = zeros(n);
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
F32 = zeros(n); F33 = zeros(n); F34 = zeros(n);
for k=1:n
    ek = E(:,k);
    F32(:,k) = f3(ek,x,x);
    F33(:,k) = f3(x,ek,x);
    F34(:,k) = f3(x,x,ek);
end
end