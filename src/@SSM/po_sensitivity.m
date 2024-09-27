function varargout = po_sensitivity(obj,W,R,DW,DR,optdof,Q,epsilon)
% po_SENSITIVITY This function calculates the sensitivity of periodic orbit. 
% In particular, it returns the sensitivity of the response amplitude at a
% given excitation frequency Omega, and also the peak position on the
% forced response curve (FRC) in reduced coordinates
%
% W:  a structure array with the SSM parameterization
% R:  a structure array with the reduced dynamics
% DW: a structure array with the sensitivity of SSM parameterization
% DR: a structure array with the sensitivity of the reduced dynamics
% optdof: dofs of response amplitude
% Q: symmetric, weight matrix
%

% setup
Lambda_E = R.Lambda_E;
gamma    = R.r21;
ftilde   = R.ftilde;
Dlamd    = DR.Dlamd;
Dgamma   = DR.Dr21;
Dftilde  = DR.Dftilde;
Omega    = obj.System.Omega;

% amplitude of periodic orbit at a given Omega
% solve for rho and theta
rhosol = rho_roots(Lambda_E,gamma,epsilon,ftilde,Omega);
thsol  = cal_theta(Lambda_E,gamma,ftilde,Omega,rhosol);

% compute periodic orbit and its amplitude
[amp,zt] = cal_pos(W,rhosol,thsol,Q,epsilon,optdof);

% solve for rho' and theta'
[Drho,Dth] = fixed_points_derivative(Lambda_E,gamma,ftilde,Dlamd,...
    Dgamma,Dftilde,rhosol,thsol,Omega,epsilon);

% compute the derivative of amplitude
Damp = cal_amp_derivative(W,DW,rhosol,thsol,Drho,Dth,Q,epsilon,optdof,amp);

% FRC peak
[rhomax,Ommax,Drhomax,DOmmax] = FRC_peak(Lambda_E,gamma,epsilon,ftilde,...
    Dlamd,Dgamma,Dftilde);

% record output
FRC = struct();
FRC.rho = rhosol;
FRC.th  = thsol;
FRC.amp = amp;
FRC.zt  = zt;
FRC.rhomax = rhomax; FRC.Ommax = Ommax;
DFRC = struct();
DFRC.Drho = Drho;
DFRC.Dth  = Dth;
DFRC.Damp = Damp;
DFRC.Drhomax = Drhomax;
DFRC.DOmmax  = DOmmax;
varargout{1} = FRC;
varargout{2} = DFRC;

end


function y = rho_roots(Lambda_E,gamma,epsilon,ftilde,x)
% RHO_ROOTS: find rho via the roots of polynomial

% setup
Rlamd = real(Lambda_E);
Ilamd = imag(Lambda_E);
Rgamma = real(gamma);
Igamma = imag(gamma);
% compute coefficients
a0 = -epsilon^2*(abs(ftilde))^2;              % rho^0
a2 = Rlamd^2+(Ilamd-x)^2;                     % rho^2
a4 = 2*Rlamd*Rgamma+2*(Ilamd-x)*Igamma;       % rho^4
a6 = Rgamma^2+Igamma^2;                       % rho^6
% find roots
coeffs = [a6 0 a4 0 a2 0 a0];
rho = roots(coeffs);
rho = real(rho(imag(rho)==0));
y   = rho(rho>0); y = y(:).';
end

function y = cal_theta(Lambda_E,gamma,ftilde,Omega,rho)
% setup
Rlamd = real(Lambda_E);
Ilamd = imag(Lambda_E);
Rgamma = real(gamma);
Igamma = imag(gamma);
Rftilde = real(ftilde);
Iftilde = imag(ftilde);
% nummerator of cos and sin
costh = -Rftilde*(Rlamd*rho+Rgamma*rho.^3)-Iftilde*((Ilamd-Omega)*rho+Igamma*rho.^3);
sinth = -Iftilde*(Rlamd*rho+Rgamma*rho.^3)+Rftilde*((Ilamd-Omega)*rho+Igamma*rho.^3);
% calculate angle
y = atan2(sinth,costh);
end

function [y,zt] = cal_pos(W,rho,th,Q,epsilon,optdof)
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
nsol = numel(rho);
zt = cell(nsol,1); amp = zeros(nsol,1);
for k=1:nsol
    rhok = rho(k); thk = th(k);
    tmp1 = ((W10*rhok+W21*rhok^3)*exp(1i*thk)+epsilon*x0)*exp(1i*tps);
    tmp2 = W20*rhok^2*exp(1i*(2*thk+2*tps))+W30*rhok^3*exp(1i*(3*thk+3*tps));
    zk   = 2*real(tmp1)+2*real(tmp2)+W11*rhok^2;
    ampk = trapz(tps,sum(zk.*(Q*zk),1))/(2*pi);
    zt{k} = zk; amp(k) = sqrt(ampk);
end
% calculate amplitude of periodic orbit
a2 = 2*W10'*Q*W10;
a4 = 2*W10'*Q*W21+2*W21'*Q*W10+2*W20'*Q*W20+W11.'*Q*W11;
a6 = 2*(W21'*Q*W21+W30'*Q*W30);
AautoSquare = a2*rho.^2+a4*rho.^4+a6*rho.^6;
AcoupSquare = 4*real(exp(1i*th).*(x0'*Q*(W10*rho+W21*rho.^3))); % rho should be in a row
AnonautoSquare = 2*x0'*Q*x0;
y = sqrt(AautoSquare+epsilon*AcoupSquare+epsilon^2*AnonautoSquare);
end

function [Drho,Dth] = fixed_points_derivative(Lambda_E,gamma,ftilde,Dlamd,Dgamma,Dftilde,rho,th,Omega,epsilon)
% setup
nsol = numel(rho);
npar = numel(Dlamd);
Drho = zeros(nsol,npar);
Dth  = zeros(nsol,npar);
Rlamd = real(Lambda_E);
Ilamd = imag(Lambda_E);
Rgamma   = real(gamma);
Igamma   = imag(gamma);
Rftilde  = real(ftilde);
Iftilde  = imag(ftilde);
Dlamd    = Dlamd(:).';
Dgamma   = Dgamma(:).';
Dftilde  = Dftilde(:).';
RDlamd   = real(Dlamd);
IDlamd   = imag(Dlamd);
RDgamma  = real(Dgamma);
IDgamma  = imag(Dgamma);
RDftilde = real(Dftilde);
IDftilde = imag(Dftilde);
% coefficient matrix
a11 = Rlamd+3*rho.^2*Rgamma;
a12 = epsilon*(-Rftilde*sin(th)+Iftilde*cos(th));
a21 = Ilamd-Omega+3*rho.^2*Igamma;
a22 = -epsilon*(Iftilde*sin(th)+Rftilde*cos(th));
% matrix inverse to compute derivatives
for k=1:nsol
    rhok = rho(k);
    thk  = th(k);
    Ak   = [a11(k) a12(k); a21(k) a22(k)];
    brho = RDlamd*rhok+RDgamma*rhok^3+epsilon*(RDftilde*cos(thk)+IDftilde*sin(thk));
    brho = -brho;
    bth  = IDlamd*rhok+IDgamma*rhok^3+epsilon*(IDftilde*cos(thk)-RDftilde*sin(thk));
    bth  = -bth;
    tmp  = Ak\[brho; bth];
    Drho(k,:) = tmp(1,:);
    Dth(k,:)  = tmp(2,:);
end
end

function y = cal_amp_derivative(W,DW,rho,th,Drho,Dth,Q,epsilon,optdof,amp)
% setup 
W10 = W.W10(optdof,:);
W20 = W.W20(optdof,:);
W11 = W.W11(optdof,:);
W30 = W.W30(optdof,:);
W21 = W.W21(optdof,:);
x0  = W.x0(optdof,:);
% derivatives
DW10 = DW.DW10(optdof,:);
DW20 = DW.DW20(optdof,:);
DW11 = DW.DW11(optdof,:);
DW30 = DW.DW30(optdof,:);
DW21 = DW.DW21(optdof,:);
Dx0  = DW.Dx0(optdof,:);

% calculate amplitude of periodic orbit
a2 = 2*W10'*Q*W10;
a4 = 2*W10'*Q*W21+2*W21'*Q*W10+2*W20'*Q*W20+W11.'*Q*W11;
a6 = 2*(W21'*Q*W21+W30'*Q*W30);
Da2 = 4*real(W10'*Q*DW10);
Da4 = 2*real(2*W10'*Q*DW21+2*W21'*Q*DW10+2*W20'*Q*DW20+W11.'*Q*DW11);
Da6 = 4*real(W21'*Q*DW21+W30'*Q*DW30);
nsol = numel(amp);
npar = size(Drho,2);
y = zeros(nsol,npar);
AnonautoDAnonauto = 2*real(x0'*Q*Dx0);
for k=1:nsol
    rhok = rho(k); thk = th(k); Drhok = Drho(k,:); Dthk = Dth(k,:);
    AautoDAauto = Da2*rhok^2+Da4*rhok^4+Da6*rhok^6+2*a2*rhok*Drhok+...
        4*a4*rhok^3*Drhok+6*a6*rhok^5*Drhok;
    AautoDAauto = AautoDAauto/2;
    tmp1 = 1i*Dthk*exp(1i*thk)*(x0'*Q*(W10*rhok+W21*rhok^3));
    tmp2 = exp(1i*thk)*(Dx0'*Q*(W10*rhok+W21*rhok^3));
    tmp3 = exp(1i*thk)*(x0'*Q*(DW10*rhok+W10*Drhok+DW21*rhok^3+3*W21*rhok^2*Drhok));
    AcoupDAcoup = 2*real(tmp1+tmp2.'+tmp3);
    y(k,:) = (AautoDAauto+epsilon*AcoupDAcoup+epsilon^2*AnonautoDAnonauto)/amp(k);
end
end

function [rhomax,Ommax,Drhomax,DOmmax] = FRC_peak(Lambda_E,gamma,epsilon,ftilde,...
    Dlamd,Dgamma,Dftilde)
% setup
Rlamd = real(Lambda_E);
Ilamd = imag(Lambda_E);
Rgamma = real(gamma);
Igamma = imag(gamma);
RDlamd = real(Dlamd);
IDlamd = imag(Dlamd);
RDgamma = real(Dgamma);
IDgamma = imag(Dgamma);

if Rlamd<0
    epsilon = -epsilon; % to ensure positive rho
end
rhomax = epsilon*abs(ftilde)/Rlamd; % linear guess
rhos   = roots([Rgamma,0,Rlamd,-epsilon*abs(ftilde)]);
[~,idx] = min(abs(rhos-rhomax));
rhomax = abs(rhos(idx));
Ommax  = Ilamd+Igamma*rhomax^2;

% sensitivity of the peak
dabsftilde = real(conj(ftilde)*Dftilde)/abs(ftilde);
Drhomax = (epsilon*dabsftilde-RDlamd*rhomax-RDgamma*rhomax^3)/(Rlamd+3*Rgamma*rhomax^2);
DOmmax  = IDlamd+IDgamma*rhomax^2+2*Igamma*rhomax*Drhomax;
end

