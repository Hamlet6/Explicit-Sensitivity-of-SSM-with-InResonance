function [W_0,R_0,W_1,R_1] = ssmtool_comp(m,c1,c2,b1,b2,f1,f2,Om)


[mass,damp,stiff,fnl,fext] = build_model(m,c1,c2,b1,b2,f1,f2);
DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl',fnl);
set(DS.Options,'Emax',2,'Nmax',4,'notation','tensor') % second-order
% Forcing
epsilon = 1e-2;
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);
% Linear Modal analysis
[V,D,W] = DS.linear_spectral_analysis();

%% SSM computation
S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','tensor');
resonant_modes = [1 2 3 4]; % choose master spectral subspace
S.choose_E(resonant_modes);
order = 3;                  % SSM expansion order
[W_0,R_0] = S.compute_whisker(order);
[W_1,R_1] = S.compute_perturbed_whisker(0,W_0,R_0,Om);
end