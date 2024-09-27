function [DW,DR]= explicit_senstivity_whisker(obj,W,R)
% Explicit senstivity of SSM upto cubic order. The current implementation
% is restricted to second order system with symmetric mass and stiffness
% matrices along with propotional damping. 

%% setup 
DW = struct();
DR = struct();
W20 = W.W20;
W11 = W.W11;
W30 = W.W30;
W21 = W.W21;
r21 = R.r21;
% check system setup
is2ndorder = obj.System.order==2; % second order system
% ispropDamp = ~obj.System.nl_damp && obj.System.Options.RayleighDamping; % damping type
ispropDamp = obj.System.Options.RayleighDamping; % damping type
assert(is2ndorder, "current implementation assumes 2nd system");
% spectrum
Lambda_E = obj.E.spectrum;
assert(numel(Lambda_E)==2,'current implementation only supports two-dimensional SSMs');
idx = 1;
if imag(Lambda_E(1))<0; idx=2; end
Lambda_E = Lambda_E(idx); midx = obj.System.n;
vE  = obj.E.basis(:,idx); vEbar = conj(vE);
uE  = obj.E.adjointBasis(:,idx);
phi = vE(1:midx);         psi = conj(uE(1:obj.System.n));

assert(abs(R.Lambda_E-Lambda_E)<1e-3*abs(R.Lambda_E),'Lambda_E is not matched');
% second-order system
M   = obj.System.M;
K   = obj.System.K;
C   = obj.System.C;
if ispropDamp
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
    kapa = -1i/(2*om*sqrt(1-zet^2));
    psi  = kapa*phi;
    KMmat = (K-lamdhat*M+lamdhat*M*phi*(phi'*M));
end
f2  = obj.System.fnl_semi{1};
f2  = @(x1,x2) f2({x1;x2});
f3  = obj.System.fnl_semi{2};
f3  = @(x1,x2,x3) f3({x1;x2;x3});

% derivative w.r.t parameters
DM  = obj.System.DM;
DK  = obj.System.DK;
nps = numel(DM); % number of parameters
if ispropDamp
    DC = cell(nps,1);
else
    DC = obj.System.DC;
end
Dfnl2 = obj.System.Dfnl2_semi;
Dfnl3 = obj.System.Dfnl3_semi;
Df2 = cell(nps,1);
Df3 = cell(nps,1);
for k = 1:nps
    Df2{k} = @(x1,x2) Dfnl2{k}({x1;x2});
    Df3{k} = @(x1,x2,x3) Dfnl3{k}({x1;x2;x3});
end

%% linear terms
Dlamd = zeros(nps,1);
Dom   = zeros(nps,1);
Dpsi  = zeros(midx,nps);
DW10  = zeros(2*midx,nps);
F11   = Lambda_E^2*M+Lambda_E*C+K; 
F12   = (2*Lambda_E*M+C)*phi;
F21   = psi.'*(2*Lambda_E*M+C);
DFmat = [F11 F12; F21 0];
for k=1:nps
    if ispropDamp
        dom   = (phi'*(DK{k}-lamdhat*DM{k})*phi)/(2*om); 
        rhs   = 2*om*dom*M*phi+lamdhat*DM{k}*phi-DK{k}*phi-0.5*lamdhat*M*phi*(phi'*DM{k}*phi);
        dphi  = KMmat\rhs;
        dzeta = (dal+dbe*om^2+2*be*om-2*zet)*dom/(2*om);
        dlamd = -(dzeta*om*Lambda_E+(zet*Lambda_E+om)*dom)/(Lambda_E+zet*om);
        dkapa = 1i/(2*om^2*(1-zet^2))*(dom*sqrt(1-zet^2)-dzeta*zet*om/sqrt(1-zet^2));
        dpsi  = dkapa*phi+kapa*dphi;
        Dom(k)   = dom;
        DC{k} = al*DM{k}+be*DK{k}+(dal*M+dbe*K)*dom;
    else
        tmp1 = Lambda_E^2*DM{k}+Lambda_E*DC{k}+DK{k};
        rhs1 = [-tmp1*phi; 0];
        sol1 = DFmat\rhs1;
        dphi = sol1(1:end-1);
        dlamd = sol1(end);
        tmp2 = (2*Lambda_E*M+C).'*psi;
        G = -tmp2*dlamd-tmp1.'*psi;
        r = -tmp2.'*dphi-2*psi.'*M*phi*dlamd-psi.'*(2*Lambda_E*DM{k}+DC{k})*phi;
        sol2 = transpose(DFmat)\[G;r];
        dpsi = sol2(1:end-1);
        assert(abs(sol2(end))<1e-4,'solution for xi is non-zero');
    end
    Dlamd(k) = dlamd;
    Dpsi(:,k)  = dpsi;
    DW10(1:midx,k) = dphi;
    DW10(midx+1:end,k) = Lambda_E*dphi+dlamd*phi;
end
obj.System.DC = DC;

%% quadratic order terms
solver = obj.Options.solver;
DW20 = zeros(2*midx,nps);
DW11 = zeros(2*midx,nps);
A20  = 4*Lambda_E^2*M+2*Lambda_E*C+K;
A11  = 4*(real(Lambda_E))^2*M+2*real(Lambda_E)*C+K;
W201 = W20(1:midx);
W111 = W11(1:midx);
phibar    = conj(phi);
RLambda_E = real(Lambda_E);
RDlamd    = real(Dlamd);
for k=1:nps
    % W20
    if obj.System.nl_damp
        dvE  = DW10(:,k);
        df2  = Df2{k}(vE,vE)+f2(dvE,vE)+f2(vE,dvE);
        dvEbar = conj(dvE);
        df2hat = Df2{k}(vE,vEbar)+f2(dvE,vEbar)+f2(vE,dvEbar)+...
            Df2{k}(vEbar,vE)+f2(dvEbar,vE)+f2(vEbar,dvE);
        df2hat = 0.5*df2hat;
    else
        dphi = DW10(1:midx,k);
        df2  = Df2{k}(phi,phi)+f2(dphi,phi)+f2(phi,dphi);
        if ispropDamp
            df2hat = df2;
        else
            dphibar = conj(dphi);
            df2hat = Df2{k}(phi,phibar)+f2(dphi,phibar)+f2(phi,dphibar)+...
                Df2{k}(phibar,phi)+f2(dphibar,phi)+f2(phibar,dphi);
            df2hat = 0.5*df2hat;            
        end
    end
    DA20 = 4*Lambda_E^2*DM{k}+2*Lambda_E*DC{k}+DK{k}+8*Lambda_E*Dlamd(k)*M+2*Dlamd(k)*C;
    DW20(1:midx,k)     = -solveinveq(A20,df2+DA20*W201,solver);
    DW20(midx+1:end,k) = 2*Dlamd(k)*W201+2*Lambda_E*DW20(1:midx,k);
    % W11
    DA11 = 4*RLambda_E^2*DM{k}+2*RLambda_E*DC{k}+DK{k}+8*RLambda_E*RDlamd(k)*M+2*RDlamd(k)*C;
    DW11(1:midx,k)     = -solveinveq(A11,2*df2hat+DA11*W111,solver);
    DW11(midx+1:end,k) = 2*RDlamd(k)*W111+2*RLambda_E*DW11(1:midx,k);    
end


%% cubic terms
DW30 = zeros(2*midx,nps);
DW21 = zeros(2*midx,nps);
A30  = 9*Lambda_E^2*M+3*Lambda_E*C+K;
lamdt = 2*Lambda_E+conj(Lambda_E);
A21  = lamdt^2*M+lamdt*C+K;
W301 = W30(1:midx);
W211 = W21(1:midx);
W212 = W21(midx+1:end);
Dr21 = zeros(nps,1);
if obj.System.nl_damp
    f21  = f2(vEbar,W20)+f2(W20,vEbar)+f2(vE,W11)+f2(W11,vE)+...
        f3(vEbar,vE,vE)+f3(vE,vEbar,vE)+f3(vE,vE,vEbar);
else
    f21  = f2(phi,W201)+f2(W201,phi)+f2(phi,W111)+f2(W111,phi)+3*f3(phi,phi,phi);
end
varho = r21/(2*RLambda_E);
for k=1:nps
    dphi = DW10(1:midx,k);
    % W30
    if obj.System.nl_damp
        dvE   = DW10(:,k);
        DW20k = DW20(:,k);
        df30  = Df2{k}(vE,W20)+f2(dvE,W20)+f2(vE,DW20k)+...
               Df2{k}(W20,vE)+f2(DW20k,vE)+f2(W20,dvE)+...
               Df3{k}(vE,vE,vE)+f3(dvE,vE,vE)+f3(vE,dvE,vE)+f3(vE,vE,dvE);
    else
        DW201 = DW20(1:midx,k);
        df30 = Df2{k}(phi,W201)+f2(dphi,W201)+f2(phi,DW201)+...
               Df2{k}(W201,phi)+f2(DW201,phi)+f2(W201,dphi)+...
               Df3{k}(phi,phi,phi)+f3(dphi,phi,phi)+f3(phi,dphi,phi)+f3(phi,phi,dphi);
    end
    DA30 = 9*Lambda_E^2*DM{k}+3*Lambda_E*DC{k}+DK{k}+18*Lambda_E*Dlamd(k)*M+3*Dlamd(k)*C;
    DW30(1:midx,k)     = -solveinveq(A30,df30+DA30*W301,solver);
    DW30(midx+1:end,k) = 3*Dlamd(k)*W301+3*Lambda_E*DW30(1:midx,k);
    % R21
    if obj.System.nl_damp
        dvEbar = conj(dvE);
        DW11k  = DW11(:,k);
        df21  = Df2{k}(vEbar,W20)+f2(dvEbar,W20)+f2(vEbar,DW20k)+...
                Df2{k}(W20,vEbar)+f2(DW20k,vEbar)+f2(W20,dvEbar)+...
                Df2{k}(vE,W11)+f2(dvE,W11)+f2(vE,DW11k)+...
                Df2{k}(W11,vE)+f2(DW11k,vE)+f2(W11,dvE)+...
                Df3{k}(vEbar,vE,vE)+f3(dvEbar,vE,vE)+f3(vEbar,dvE,vE)+f3(vEbar,vE,dvE)+...
                Df3{k}(vE,vEbar,vE)+f3(dvE,vEbar,vE)+f3(vE,dvEbar,vE)+f3(vE,vEbar,dvE)+...
                Df3{k}(vE,vE,vEbar)+f3(dvE,vE,vEbar)+f3(vE,dvE,vEbar)+f3(vE,vE,dvEbar);
    else
        DW111 = DW11(1:midx,k);
        if ispropDamp
            df21  = Df2{k}(phi,W201)+f2(dphi,W201)+f2(phi,DW201)+...
                    Df2{k}(W201,phi)+f2(DW201,phi)+f2(W201,dphi)+...
                    Df2{k}(phi,W111)+f2(dphi,W111)+f2(phi,DW111)+...
                    Df2{k}(W111,phi)+f2(DW111,phi)+f2(W111,dphi)+...
                    3*(Df3{k}(phi,phi,phi)+f3(dphi,phi,phi)+f3(phi,dphi,phi)+f3(phi,phi,dphi));
        else
            dphibar = conj(dphi);
            df21  = Df2{k}(phibar,W201)+f2(dphibar,W201)+f2(phibar,DW201)+...
                    Df2{k}(W201,phibar)+f2(DW201,phibar)+f2(W201,dphibar)+...
                    Df2{k}(phi,W111)+f2(dphi,W111)+f2(phi,DW111)+...
                    Df2{k}(W111,phi)+f2(DW111,phi)+f2(W111,dphi)+...
                    Df3{k}(phibar,phi,phi)+f3(dphibar,phi,phi)+f3(phibar,dphi,phi)+f3(phibar,phi,dphi)+...
                    Df3{k}(phi,phibar,phi)+f3(dphi,phibar,phi)+f3(phi,dphibar,phi)+f3(phi,phibar,dphi)+...
                    Df3{k}(phi,phi,phibar)+f3(dphi,phi,phibar)+f3(phi,dphi,phibar)+f3(phi,phi,dphibar);
        end
    end
    Dr21(k) = -Dpsi(:,k).'*f21-psi.'*df21;
    % W21
    dlamdt = 2*Dlamd(k)+conj(Dlamd(k));
    DA21   = lamdt^2*DM{k}+lamdt*DC{k}+DK{k}+2*lamdt*dlamdt*M+dlamdt*C;
    dvarho = Dr21(k)/(2*RLambda_E)-r21*RDlamd(k)/(2*RLambda_E^2);
    DW21(1:midx,k)     = -solveinveq(A21,df21+DA21*(W211+varho*phi),solver)-dvarho*phi-varho*dphi;
    tmp1               = dlamdt*f21+lamdt*df21+DA21*(W212+varho*Lambda_E*phi);
    tmp2               = dvarho*Lambda_E*phi+varho*Dlamd(k)*phi+varho*Lambda_E*dphi;
    DW21(midx+1:end,k) = -solveinveq(A21,tmp1,solver)-tmp2;
end

DR.Dom = Dom;
DR.Dpsi = Dpsi;

% leading-order of non-autonomous part
fext = obj.System.fext;
Dftilde = zeros(nps,1);
Dx0     = zeros(2*midx,nps);
if ~isempty(fext)
    ncol = size(fext,2);
    if ncol>1
        assert(ncol==2,'only a pair of forcing is allowed');
        fext = 2*fext(:,1);
    end
    dfext = obj.System.dfext;
    Om  = obj.System.Omega;
    AOm = [-K-1j*Om*C -1j*Om*M; -1j*Om*M M];
    ftilde = 0.5*psi.'*fext;
    bOm = [(C+Lambda_E*M)*phi*ftilde-0.5*fext; M*phi*ftilde];
    x0  = solveinveq(AOm,bOm,solver);   
    for k=1:nps
        dphi = DW10(1:midx,k);
        dpsi = Dpsi(:,k);
        DAOm = [-DK{k}-1j*Om*DC{k} -1j*Om*DM{k}; -1j*Om*DM{k} DM{k}];
        Dftilde(k) = 0.5*(dpsi.'*fext+psi.'*dfext{k});
        tmp1 = (DC{k}+Dlamd(k)*M+Lambda_E*DM{k})*phi*ftilde+...
            (C+Lambda_E*M)*(dphi*ftilde+phi*Dftilde(k))-0.5*dfext{k};
        tmp2 = DM{k}*phi*ftilde+M*dphi*ftilde+M*phi*Dftilde(k);
        DbOm = [tmp1; tmp2];
        Dx0(:,k) = solveinveq(AOm, DbOm-DAOm*x0,solver);
    end
    DR.Dftilde = Dftilde;
    DW.Dx0 = Dx0;
end

DW.DW10 = DW10;
DW.DW20 = DW20;
DW.DW11 = DW11;
DW.DW30 = DW30;
DW.DW21 = DW21;
DR.Dlamd = Dlamd;
DR.Dr21 = Dr21;
end