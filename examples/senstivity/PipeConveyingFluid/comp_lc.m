function pLC = comp_lc(n,flowspeed,beta,alpha,lo)

% create model and setup
[M,C,K,fnl,~,DM,DC,DK,dfnl2,dfnl3,~] = build_model_semi(n,flowspeed,beta,alpha,'clamped-free');
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl_semi',fnl,'nl_damp',true);
set(DS.Options,'Emax',10,'Nmax',10,'notation','multiindex');
set(DS.Options,'RayleighDamping',false);

[V,D,W] = DS.linear_spectral_analysis();
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
pLC = S.LC_sensitivity(We,Re,DW,DR,1:2*n,1,0);

end