function [W_0,R_0,W_1,R_1] = ssmtool_comp(mu,f,n,Om)


[mass,damp,stiff,fnl,fext] = build_model(mu,f,n);
DS = DynamicalSystem();
% set(DS,'M',mass,'C',damp,'K',stiff,'fnl',fnl);
B = [damp,mass;mass,zeros(n,n)];
A = [-stiff,zeros(n,n);zeros(n,n),mass];
% F = [-fnl;zeros(n,1)];
set(DS,'B',B,'A',A,'fnl',fnl);
set(DS.Options,'Emax',4,'Nmax',10,'notation','tensor') % first-order
% Forcing
epsilon = 5e-2;
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