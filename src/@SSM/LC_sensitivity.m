function varargout = LC_sensitivity(obj,W,R,DW,DR,outdof,pidx,pfracs)
% FRC_SENSITIVITY This function calculates the sensitivity of limit cycles
% (LC). In particular, it returns the perturbed LCs when
% the pidx-th parameter is perturbed with pfractions. Here we assume pidx
% is a scalar, yet the pfractions can be a vector.
%
% W:  a structure array with the SSM parameterization
% R:  a structure array with the reduced dynamics
% DW: a structure array with the sensitivity of SSM parameterization
% DR: a structure array with the sensitivity of the reduced dynamics
% outdof:  dofs of response amplitude
% pidx:    index for perturbed parameters
% pfracs:  perturbation fractions (should be much less than one)

% setup
Lambda_E = R.Lambda_E;
gamma    = R.r21;
Dlamd    = DR.Dlamd;
Dgamma   = DR.Dr21;

nfrac = numel(pfracs);
pLCs  = cell(nfrac,1);
for k=1:nfrac
    % compute perturbed coefficients - reduced dynamics
    lamdk   = Lambda_E+Dlamd(pidx)*pfracs(k);
    gammak  = gamma+Dgamma(pidx)*pfracs(k);
    Rlamd   = real(lamdk);
    Ilamd   = imag(lamdk); 
    Rgamma  = real(gammak);
    Igamma  = imag(gammak);
    assert(Rlamd*Rgamma<0,'limit cycle does not exist');
    % LC in reduced coordinates
    rhoast = sqrt(-Rlamd/Rgamma);
    omast  = Ilamd+Igamma*rhoast^2;
    % compute perturbed coefficients - expansion coefficients
    Wk = struct();
    Wk.W10 = W.W10+DW.DW10(:,pidx)*pfracs(k);
    Wk.W11 = W.W11+DW.DW11(:,pidx)*pfracs(k);
    Wk.W20 = W.W20+DW.DW20(:,pidx)*pfracs(k);
    Wk.W30 = W.W30+DW.DW30(:,pidx)*pfracs(k);
    Wk.W21 = W.W21+DW.DW21(:,pidx)*pfracs(k);
    % BB in physical coordinates
    zoutk = cal_pos(Wk,rhoast,outdof);
    LCk = struct();
    LCk.rhoast = rhoast;
    LCk.omast  = omast;
    LCk.zout   = zoutk;
    pLCs{k}  = LCk;
end

% compute periodic orbit and its amplitude
varargout{1} = pLCs;
end


function zout = cal_pos(W,rho,outdof)
% setup
W10 = W.W10(outdof,:);
W20 = W.W20(outdof,:);
W11 = W.W11(outdof,:);
W30 = W.W30(outdof,:);
W21 = W.W21(outdof,:);
% calculate periodic orbit
nps = 1280;
tps = linspace(0,2*pi,nps);
tmp1 = (W10*rho+W21*rho^3)*exp(1i*tps);
tmp2 = W20*rho^2*exp(1i*2*tps)+W30*rho^3*exp(1i*3*tps);
zout  = 2*real(tmp1)+2*real(tmp2)+W11*rho^2;
end
