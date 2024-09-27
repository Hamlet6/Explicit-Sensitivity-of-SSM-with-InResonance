function ssm = comp_ssm(n,alpha,gamma,k1,kf,lo)

% create model and setup
[mass,damp,gyro,stiff,fnl,~,~,~,~,~] = build_model_semi(n,alpha,gamma,k1,kf);
DS = DynamicalSystem();
set(DS,'M',mass,'C',damp+gyro,'K',stiff,'fnl_semi',fnl,'nl_damp',true);
set(DS.Options,'Emax',6,'Nmax',10,'notation','multiindex');
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
% SSM computation
S = SSM(DS);
S.choose_E(resonant_modes);
[W10,W20,W11,W30,W21,r21] = S.explicit_whisker();

ssm = struct();
ssm.W10 = W10;
ssm.W20 = W20;
ssm.W11 = W11;
ssm.W21 = W21;
ssm.W30 = W30;
ssm.r21 = r21;
ssm.lamd = S.System.spectrum.Lambda(1);
ssm.psi  = conj(uE(1:n));

end


