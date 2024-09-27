function varargout = FRC_sensitivity(obj,W,R,DW,DR,optdof,Q,epsilon,freqRange,pidx,pfracs,varargin)
% FRC_SENSITIVITY This function calculates the sensitivity of forced
% response curve (FRC). In particular, it returns the perturbed FRCs when
% the pidx-th parameter is perturbed with pfractions. Here we assume pidx
% is a scalar, yet the pfractions can be a vector.
%
% W:  a structure array with the SSM parameterization
% R:  a structure array with the reduced dynamics
% DW: a structure array with the sensitivity of SSM parameterization
% DR: a structure array with the sensitivity of the reduced dynamics
% optdof:  dofs of response amplitude
% Q:       symmetric, weight matrix
% epsilon: forcing amplitude
% pidx:    index for perturbed parameters
% pfracs:  perturbation fractions (should be much less than one)

% setup
Lambda_E = R.Lambda_E;
gamma    = R.r21;
ftilde   = R.ftilde;
Dlamd    = DR.Dlamd;
Dgamma   = DR.Dr21;
Dftilde  = DR.Dftilde;

nfrac = numel(pfracs);
ndof  = numel(optdof);
pFRCs = cell(nfrac,1);

for k=1:nfrac
    % compute perturbed coefficients - reduced dynamics
    lamdk   = Lambda_E+Dlamd(pidx)*pfracs(k);
    gammak  = gamma+Dgamma(pidx)*pfracs(k);
    ftildek = ftilde+Dftilde(pidx)*pfracs(k);
    % FRC in (rho, omega)
    FRCk = cubic_FRC_rho_om(lamdk,gammak,ftildek,epsilon,freqRange,varargin{:});
    % compute perturbed coefficients - expansion coefficients
    Wk = struct();
    Wk.W10 = W.W10+DW.DW10(:,pidx)*pfracs(k);
    Wk.W11 = W.W11+DW.DW11(:,pidx)*pfracs(k);
    Wk.W20 = W.W20+DW.DW20(:,pidx)*pfracs(k);
    Wk.W30 = W.W30+DW.DW30(:,pidx)*pfracs(k);
    Wk.W21 = W.W21+DW.DW21(:,pidx)*pfracs(k);
    if obj.Options.contribNonAuto
        switch obj.FRCOptions.interpMethod
            case 'constant'
                % constant interpolation via the middle point
                om1 = (freqRange(1)+freqRange(2))/2;
                [x01,Dx01] = cal_leading_nonautonmous_explict(obj,DW,DR,om1,pidx);
                x01p   = x01+Dx01*pfracs(k);
                x0func = @(x) x01p;
            case 'linear'
                % linear interpolation
                om1 = freqRange(1);
                om2 = freqRange(2);
                [x01,Dx01] = cal_leading_nonautonmous_explict(obj,DW,DR,om1,pidx);
                [x02,Dx02] = cal_leading_nonautonmous_explict(obj,DW,DR,om2,pidx);
                x01p = x01+Dx01*pfracs(k);
                x02p = x02+Dx02*pfracs(k);
                x0func = @(x) x01p*(x-om2)/(om1-om2)+x02p*(x-om1)/(om2-om1);                
            case 'quadratic'
                % quadratic interpolation via three points
                om1 = freqRange(1);
                om2 = (freqRange(1)+freqRange(2))/2;
                om3 = freqRange(2);
                [x01,Dx01] = cal_leading_nonautonmous_explict(obj,DW,DR,om1,pidx);
                [x02,Dx02] = cal_leading_nonautonmous_explict(obj,DW,DR,om2,pidx);
                [x03,Dx03] = cal_leading_nonautonmous_explict(obj,DW,DR,om3,pidx);
                x01p = x01+Dx01*pfracs(k);
                x02p = x02+Dx02*pfracs(k);
                x03p = x03+Dx03*pfracs(k);
                x0func = @(x) [x01p x02p x03p]*Lquadfunc(x,om1,om2,om3);
        end
    else
        x0func = @(x) zeros(size(Wk.W10));
    end
    % FRC in physical coordinates
    if ~iscell(FRCk) % without isola
        npts = numel(FRCk.rho);
        ampL2  = zeros(npts,1); zt = cell(npts,1); ampLinf = zeros(ndof,npts);
        for j=1:npts
            Wk.x0 = x0func(FRCk.om(j));
            [ampL2(j),zt{j},ampLinf(:,j)] = cal_pos(Wk,FRCk.rho(j),FRCk.th(j),Q,epsilon,optdof);
        end
        FRCk.ampL2   = real(ampL2);
        FRCk.ampLinf = ampLinf;
        FRCk.zout = zt;
    else            % with isola
        for i=1:2
            npts = numel(FRCk{i}.rho);
            ampL2  = zeros(npts,1); zt = cell(npts,1); ampLinf = zeros(ndof,npts);
            for j=1:npts
                Wk.x0 = x0func(FRCk{i}.om(j));
                [ampL2(j),zt{j},ampLinf(:,j)] = cal_pos(Wk,FRCk{i}.rho(j),FRCk{i}.th(j),Q,epsilon,optdof);
            end
            FRCk{i}.ampL2   = real(ampL2);
            FRCk{i}.ampLinf = ampLinf;
            FRCk{i}.zout = zt;
        end
    end
    pFRCs{k}  = FRCk;
end

% compute periodic orbit and its amplitude
varargout{1} = pFRCs;
end


function [y,zt,varargout] = cal_pos(W,rho,th,Q,epsilon,optdof)
% setup
W10 = W.W10(optdof,:);
W20 = W.W20(optdof,:);
W11 = W.W11(optdof,:);
W30 = W.W30(optdof,:);
W21 = W.W21(optdof,:);
x0  = W.x0(optdof,:);
% calculate periodic orbit
nps = 1280;
tps = linspace(0,2*pi,nps);
nsol = numel(rho); ndof = numel(optdof);
zt = cell(nsol,1); ampL2 = zeros(nsol,1); ampLinf = zeros(ndof,nsol);
for k=1:nsol
    rhok = rho(k); thk = th(k);
    tmp1 = ((W10*rhok+W21*rhok^3)*exp(1i*thk)+epsilon*x0)*exp(1i*tps);
    tmp2 = W20*rhok^2*exp(1i*(2*thk+2*tps))+W30*rhok^3*exp(1i*(3*thk+3*tps));
    zk   = 2*real(tmp1)+2*real(tmp2)+W11*rhok^2;
    ampk = trapz(tps,sum(zk.*(Q*zk),1))/(2*pi); ampL2(k) = sqrt(ampk); % L2-norm
    ampLinf(:,k) = max(abs(zk),[],2);
    zt{k} = zk; 
end
% calculate amplitude of periodic orbit
a2 = 2*W10'*Q*W10;
a4 = 2*W10'*Q*W21+2*W21'*Q*W10+2*W20'*Q*W20+W11.'*Q*W11;
a6 = 2*(W21'*Q*W21+W30'*Q*W30);
AautoSquare = a2*rho.^2+a4*rho.^4+a6*rho.^6;
AcoupSquare = 4*real(exp(1i*th).*(x0'*Q*(W10*rho+W21*rho.^3))); % rho should be in a row
AnonautoSquare = 2*x0'*Q*x0;
y = sqrt(AautoSquare+epsilon*AcoupSquare+epsilon^2*AnonautoSquare);
varargout{1} = ampLinf; varargout{2} = ampL2;
end



function y = cubic_FRC_rho_om(lamd,gamma,ftilde,epsilon,freqRange,varargin)
% setup
Rlamd  = real(lamd);
Ilamd  = imag(lamd);
Rgamma = real(gamma);
Igamma = imag(gamma);
Rftilde = real(ftilde);
Iftilde = imag(ftilde);

% determine bounds for rho -> [rhomin1 -- rhomax -- rhomin2]
% upper bounds
rhomax = epsilon*abs(ftilde)/Rlamd;    % linear guess
rhorts = roots([Rgamma,0,Rlamd,-epsilon*abs(ftilde)]);
rhorts = abs(rhorts(imag(rhorts)==0)); % only positive roots matter
[~,idx] = min(abs(rhorts-rhomax));
rhomax  = rhorts(idx);
% lower bound 1
opts  = optimset('TolFun',1e-8,'TolX',1e-8);
famp  = epsilon*abs(ftilde);
Ommin = freqRange(1);
omfunc1 = @(r) Ilamd+Igamma*r.^2-real(sqrt(famp^2-(Rlamd*r+Rgamma*r.^3).^2)./r);
rhomin1 = famp/sqrt((Ilamd-Ommin)^2+Rlamd^2);
rhomin1 = fminsearch(@(x)abs(omfunc1(x)-Ommin), rhomin1, opts);
% lower bound2
Ommax = freqRange(2);
omfunc2 = @(r) Ilamd+Igamma*r.^2+real(sqrt(famp^2-(Rlamd*r+Rgamma*r.^3).^2)./r);
rhomin2 = famp/sqrt((Ommax-Ilamd)^2+Rlamd^2);
rhomin2 = fminsearch(@(x)abs(omfunc2(x)-Ommax), rhomin2, opts);

% compute FRC via as two segments
rhop = linspace(rhomin1,rhomax,150);
Omp  = omfunc1(rhop);
rhom = linspace(rhomax,rhomin2,150);
Omm  = omfunc2(rhom);
oms  = [Omp Omm];
rhos = [rhop rhom];

% account for isola
isola = false;
if numel(varargin)>0
    isola = strcmp(varargin{1},'isola');
    if isola && Rlamd*Rgamma<0
        isola = true; idmain = numel(rhos);
        rhoisola = rhorts(setdiff([1 2 3],idx));
        rholb  = min(rhoisola);
        rhoub  = max(rhoisola);
        rhoiso = linspace(rholb,rhoub,50);
        om1iso = omfunc1(rhoiso);
        om2iso = omfunc2(rhoiso);
        rhos = [rhos rhoiso flip(rhoiso)];
        oms  = [oms om1iso flip(om2iso)];
    end
end

% determine the stability of fixed points
npts = numel(oms);
st   = false(npts,1);
a11  = Rlamd+3*Rgamma*rhos.^2;
a12  = -(Ilamd-oms+Igamma*rhos.^2).*rhos;
a21  = 2*Igamma*rhos+(Ilamd-oms+Igamma*rhos.^2)./rhos;
a22  = Rlamd+Rgamma*rhos.^2;
traceJ = a11+a22;
detJ = a11.*a22-a12.*a21;
for k=1:npts
    st(k) = traceJ(k)<0 && detJ(k)>0;
end

% compute phase angle
costh = -Rftilde*(Rlamd*rhos+Rgamma*rhos.^3)-Iftilde*((Ilamd-oms).*rhos+Igamma*rhos.^3);
sinth = -Iftilde*(Rlamd*rhos+Rgamma*rhos.^3)+Rftilde*((Ilamd-oms).*rhos+Igamma*rhos.^3);
th = atan2(sinth,costh);

% record output
if ~isola
    y = struct();
    y.rho = rhos; y.om = oms; y.st = st; y.th = th;
else
    frcmain = struct();
    frcmain.rho = rhos(1:idmain); frcmain.om = oms(1:idmain);
    frcmain.st  = st(1:idmain);   frcmain.th = th(1:idmain);
    frcisola = struct();
    frcisola.rho = rhos(idmain+1:end); frcisola.om = oms(idmain+1:end);
    frcisola.st  = st(idmain+1:end);   frcisola.th = th(idmain+1:end);
    y = {frcmain,frcisola};
end
end


function [x0,Dx0] = cal_leading_nonautonmous_explict(obj,DW,DR,Om,pidx)
% leading-order of non-autonomous part
solver = obj.Options.solver;
Lambda_E = obj.E.spectrum;
assert(numel(Lambda_E)==2,'current implementation only supports two-dimensional SSMs');
idx = 1; if imag(Lambda_E(1))<0; idx=2; end
Lambda_E = Lambda_E(idx);
vE  = obj.E.basis(:,idx); 
uE  = obj.E.adjointBasis(:,idx);
phi = vE(1:obj.System.n); psi = conj(uE(1:obj.System.n));
M   = obj.System.M;
C   = obj.System.C;
K   = obj.System.K;

fext = obj.System.fext;
ncol = size(fext,2);
if ncol>1
    assert(ncol==2,'only a pair of forcing is allowed');
    fext = 2*fext(:,1);
end
AOm = [-K-1j*Om*C -1j*Om*M; -1j*Om*M M];
ftilde = 0.5*psi.'*fext;
bOm = [(C+Lambda_E*M)*phi*ftilde-0.5*fext; M*phi*ftilde];
x0  = solveinveq(AOm,bOm,solver);

DM  = obj.System.DM;
DC  = obj.System.DC;
DK  = obj.System.DK;
Dlamd = DR.Dlamd;
dfext = obj.System.dfext;   
dphi = DW.DW10(1:obj.System.n,pidx);
dpsi = DR.Dpsi(:,pidx);
DAOm = [-DK{pidx}-1j*Om*DC{pidx} -1j*Om*DM{pidx}; -1j*Om*DM{pidx} DM{pidx}];
Dftilde = 0.5*(dpsi.'*fext+psi.'*dfext{pidx});
tmp1 = (DC{pidx}+Dlamd(pidx)*M+Lambda_E*DM{pidx})*phi*ftilde+...
    (C+Lambda_E*M)*(dphi*ftilde+phi*Dftilde)-0.5*dfext{pidx};
tmp2 = DM{pidx}*phi*ftilde+M*dphi*ftilde+M*phi*Dftilde;
DbOm = [tmp1; tmp2];
Dx0  = solveinveq(AOm, DbOm-DAOm*x0,solver);
end

function y = Lquadfunc(x,om1,om2,om3)
% evaluation of Lagrangian basis functions at x
y    = zeros(3,1);
y(1) = (x-om2)*(x-om3)/((om1-om2)*(om1-om3));
y(2) = (x-om1)*(x-om3)/((om2-om1)*(om2-om3));
y(3) = (x-om1)*(x-om2)/((om3-om1)*(om3-om2));
end

