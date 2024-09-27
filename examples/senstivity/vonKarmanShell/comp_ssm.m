function ssm = comp_ssm(ngrids,E,rho,Om,type)

% create model and setup
[M,K,fnl,fext,~,~,~,~,~,~] = build_model_semi(ngrids,E,rho);
DS = DynamicalSystem();
switch type
    case 'type1'
        al = 0.402153037834110;
        be = 8.631167461080214e-06;
        C = al*M+be*K;
    case 'type2'
        disp('With fixed damping ratio for the first mode');
        lamd = eigs(K,M,1,'smallestabs'); 
        zeta = 0.02;
        om = sqrt(lamd);
        al = zeta*om; 
        be = zeta/om;
        C  = al*M+be*K;
end
set(DS,'M',M,'C',C,'K',K,'fnl_semi',fnl);
set(DS.Options,'Emax',3,'Nmax',6)
set(DS,'fext',fext,'Omega',Om);
[V,D,W] = DS.linear_spectral_analysis();

S = SSM(DS);
set(S.Options,'solver','backslash');
resonant_modes = [1 2]; % choose master spectral subspace
S.choose_E(resonant_modes);
[We,Re] = S.explicit_whisker();

ssm = struct();
ssm.W10 = We.W10;
ssm.W20 = We.W20;
ssm.W11 = We.W11;
ssm.W21 = We.W21;
ssm.W30 = We.W30;
ssm.psi = We.psi;
ssm.r21 = Re.r21;
ssm.lamd = S.System.spectrum.Lambda(1);
ssm.x0   = We.x0;
ssm.ftilde = Re.ftilde;
end


