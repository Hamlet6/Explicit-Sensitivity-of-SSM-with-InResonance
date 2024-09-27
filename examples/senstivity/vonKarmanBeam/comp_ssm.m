function ssm = comp_ssm(kLin,kNon)

% create model and setup
nElements = 100;
[M,~,K,fnl,~,~,~] = build_model_semi(nElements,kLin,kNon);
om = eigs(K,M,1,'smallestabs');
om = sqrt(om);
zeta = 0.01;
C    = 2*zeta*om*M;
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl_semi',fnl);
set(DS.Options,'Emax',1,'Nmax',10,'notation','multiindex');

% SSM computation
S = SSM(DS);
resonant_modes = [1 2]; % choose master spectral subspace
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

end


