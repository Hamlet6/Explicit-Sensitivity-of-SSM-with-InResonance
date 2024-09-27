function [W_0,R_0,W_1,R_1] = ssmtool_comp(nmodes,flowspeed,beta,alpha,Om)

[M, C, K, fnl, fext] = build_model(nmodes,flowspeed,beta,0,0,alpha,'clamped-free','nonlinear_damp');

% Create dynamical system 
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
epsilon = 1e-2;
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);
set(DS.Options,'Emax',10,'Nmax',10,'notation','multiindex')
set(DS.Options,'RayleighDamping',false)
[V,D,W] = DS.linear_spectral_analysis();

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