function varargout = BB_sensitivity(obj,W,R,DW,DR,outdof,freqRange,pidx,pfracs)
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
Dlamd    = DR.Dlamd;
Dgamma   = DR.Dr21;

nfrac = numel(pfracs);
if isnumeric(outdof)
    ndof = numel(outdof);
else
    ndof = numel(outdof(W.W10));
end
pBBs  = cell(nfrac,1);

for k=1:nfrac
    % compute perturbed coefficients - reduced dynamics
    lamdk   = Lambda_E+Dlamd(pidx)*pfracs(k);
    gammak  = gamma+Dgamma(pidx)*pfracs(k);
    % BB in (rho, omega)
    BBk = cubic_BB_rho_om(lamdk,gammak,freqRange);
    % compute perturbed coefficients - expansion coefficients
    Wk = struct();
    Wk.W10 = W.W10+DW.DW10(:,pidx)*pfracs(k);
    Wk.W11 = W.W11+DW.DW11(:,pidx)*pfracs(k);
    Wk.W20 = W.W20+DW.DW20(:,pidx)*pfracs(k);
    Wk.W30 = W.W30+DW.DW30(:,pidx)*pfracs(k);
    Wk.W21 = W.W21+DW.DW21(:,pidx)*pfracs(k);
    % BB in physical coordinates
    npts = numel(BBk.rho);
    ampLinf = zeros(ndof,npts);
    for j=1:npts
        ampLinf(:,j) = cal_pos(Wk,BBk.rho(j),outdof);
    end
    BBk.ampLinf = ampLinf;
    pBBs{k}  = BBk;
end

% compute periodic orbit and its amplitude
varargout{1} = pBBs;
end


function ampLinf = cal_pos(W,rho,outdof)
% setup
if isnumeric(outdof)
    uidx = outdof; ndof = numel(outdof);
else
    uidx = 1:numel(W.W10); ndof = numel(outdof(W.W10));
end
W10 = W.W10(uidx,:);
W20 = W.W20(uidx,:);
W11 = W.W11(uidx,:);
W30 = W.W30(uidx,:);
W21 = W.W21(uidx,:);
% calculate periodic orbit
nps = 1280;
tps = linspace(0,2*pi,nps);
nsol = numel(rho); 
ampLinf = zeros(ndof,nsol);
for k=1:nsol
    rhok = rho(k); 
    tmp1 = (W10*rhok+W21*rhok^3)*exp(1i*tps);
    tmp2 = W20*rhok^2*exp(1i*2*tps)+W30*rhok^3*exp(1i*3*tps);
    zk   = 2*real(tmp1)+2*real(tmp2)+W11*rhok^2;
    if isa(outdof,'function_handle')
        zk = outdof(zk);
    end
    ampLinf(:,k) = max(abs(zk),[],2);
end
end

function y = cubic_BB_rho_om(lamd,gamma,freqRange)
% setup
Ilamd  = imag(lamd);
Igamma = imag(gamma);
% determine upper bound for rho
if Igamma>0 
    rhomax = sqrt((freqRange(2)-Ilamd)/Igamma);
else
    rhomax = sqrt((freqRange(1)-Ilamd)/Igamma);
end
% samplings for rho
rhos = linspace(0,rhomax,101);
% omegas
oms  = Ilamd+Igamma*rhos.^2;
y = struct();
y.rho = rhos; y.om = oms;
end
