function [frce,varargout] = comp_frc(n,flowspeed,beta,alpha,lo,Omega,Q,epsilon,varargin)

% create model and setup
[M,C,K,fnl,fext,DM,DC,DK,dfnl2,dfnl3,dfext] = build_model_semi(n,flowspeed,beta,alpha,'clamped-free');
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl_semi',fnl,'nl_damp',true,'fext',fext);
set(DS.Options,'Emax',10,'Nmax',10,'notation','multiindex');
set(DS.Options,'RayleighDamping',false);

[V,D,W] = DS.linear_spectral_analysis();
set(DS,'dfext',dfext,'Omega',Omega);
resonant_modes = [1 2]; % choose master spectral subspace
rmode = resonant_modes(1);
% imposing normalization condition
phi = V(1:n,rmode);
phi = phi./(lo.'*phi);
vE = [phi; D(rmode)*phi];
mu = W(:,rmode)'*DS.B*vE;
uE = W(:,rmode)./(mu');

VE = [vE conj(vE)];
WE = [uE conj(uE)];

WE'*DS.B*VE

WE'*DS.A*VE-diag(D(resonant_modes))

DS.spectrum.V = VE;
DS.spectrum.W = WE;

% explict SSM
S = SSM(DS);
S.choose_E(resonant_modes);
[We,Re] = S.explicit_whisker();
% explicit derivative
set(DS,'DM',DM,'DC',DC,'DK',DK);
set(DS,'Dfnl2_semi',dfnl2,'Dfnl3_semi',dfnl3);
[DW,DR] = S.explicit_senstivity_whisker(We,Re);
optdof = 1:n; 
[frce,~] = S.po_sensitivity(We,Re,DW,DR,optdof,Q,epsilon);
if numel(varargin)>0
    freqRange = varargin{1};
    if numel(varargin)==2
        pFRC = S.FRC_sensitivity(We,Re,DW,DR,optdof,Q,epsilon,freqRange,1,0,'isola');
    else
        pFRC = S.FRC_sensitivity(We,Re,DW,DR,optdof,Q,epsilon,freqRange,1,0);
    end
    varargout{1} = pFRC;
end


end