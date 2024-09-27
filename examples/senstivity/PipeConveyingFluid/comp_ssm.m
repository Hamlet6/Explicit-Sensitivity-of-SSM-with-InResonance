function ssm = comp_ssm(n,flowspeed,beta,alpha,lo,Omega)

% create model and setup
[M,C,K,fnl,fext,~,~,~,~,~,~] = build_model_semi(n,flowspeed,beta,alpha,'clamped-free');
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl_semi',fnl,'nl_damp',true,'fext',fext);
set(DS.Options,'Emax',10,'Nmax',10,'notation','multiindex');
set(DS.Options,'RayleighDamping',false);

[V,D,W] = DS.linear_spectral_analysis();
set(DS,'fext',fext,'Omega',Omega);
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
[We,Re] = S.explicit_whisker();

ssm = struct();
ssm.W10 = We.W10;
ssm.W20 = We.W20;
ssm.W11 = We.W11;
ssm.W21 = We.W21;
ssm.W30 = We.W30;
ssm.r21 = Re.r21;
ssm.lamd = S.System.spectrum.Lambda(1);
ssm.psi  = conj(uE(1:n));
ssm.x0   = We.x0;
ssm.ftilde = Re.ftilde;
end


