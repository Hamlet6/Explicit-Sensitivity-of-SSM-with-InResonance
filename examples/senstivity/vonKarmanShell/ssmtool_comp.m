function [W_0,R_0,W_1,R_1] = ssmtool_comp(nDiscretization,E,rho,dampType,Om)

% build model
[M,C,K,fnl,fext,~] = build_model(nDiscretization,E,rho,dampType);

% Create dynamical system 
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
epsilon = 1e-2;
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);
set(DS.Options,'Emax',3,'Nmax',6,'notation','multiindex')
[V,D,W] = DS.linear_spectral_analysis();

% recale Vs and Ws such that phi'Mphi = I
phi = DS.spectrum.V(1:size(M,1),1);
tmp = sqrt(phi'*M*phi);
DS.spectrum.V = V/tmp;
DS.spectrum.W = W*tmp;

% Choose Master subspace (perform resonance analysis)
S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','multiindex');
order = 3;
resonant_modes = [1,2]; 
S.choose_E(resonant_modes);

% autonomous part
[W_0,R_0] = S.compute_whisker(order);

% non-autonomous part
[W_1,R_1] = S.compute_perturbed_whisker(0,W_0,R_0,Om);

end