function ssm = comp_ssm(m,b1,b2,zeta)

% create model and setup
[mass,damp,stiff,fnl,DM,DK,dfnl2,dfnl3] = build_model_sens_direct(m,b1,b2,zeta);
DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl_semi',fnl);
set(DS.Options,'Emax',1,'Nmax',2)
set(DS,'fext',[m;b1*b2],'Omega',0.98);
[V,D,W] = DS.linear_spectral_analysis();

% SSM computation
S = SSM(DS);
resonant_modes = [1 2]; % choose master spectral subspace
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
ssm.psi  = conj(W(1:DS.n,1));
ssm.x0   = We.x0;
ssm.ftilde = Re.ftilde;

end


