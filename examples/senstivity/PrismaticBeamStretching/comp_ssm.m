function ssm = comp_ssm(b,h,L,rho,E,P)

% create model and setup
[mass,stiff,fnl,~,~,~,~] = build_model_semi(b,h,L,rho,E,P);
om   = eigs(stiff,mass,1,'smallestabs');
om   = sqrt(om);
zeta = 0.01;
damp = 2*zeta*om*mass;
DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl_semi',fnl);
set(DS.Options,'Emax',1,'Nmax',2)

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


