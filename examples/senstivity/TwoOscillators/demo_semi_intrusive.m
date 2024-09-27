%% 
% Consider the following system with 2 dofs
% 
% $$\ddot{x}_1+c_1\dot{x}_1+x_1+b_1x_1x_2+x_1^3+x_1x_2^2=\epsilon f_1\cos\Omega 
% t\\\ddot{x}_2+c_2\dot{x}_2+2x_2+b_2x_1^2+0.3x_1x_2^2+1.6x_2^3=\epsilon f_2\cos\Omega 
% t$$

clear all
%% Setup model

m = 1;
c1 = 5e-3;
c2 = 1e-2;
b1 = 1;
b2 = 1;
f1 = 1;
f2 = 0;
[mass,damp,stiff,fnl,fext] = build_model_semi(c1,c2,b1,b2,f1,f2);

DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl_semi',fnl,'F_semi_sym',false);
set(DS.Options,'Emax',1,'Nmax',2,'notation','multiindex')
% set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(DS.Options,'Intrusion','semi')


% Forcing
epsilon = 1e-2;
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);
% Linear Modal analysis

[V,D,W] = DS.linear_spectral_analysis();
%% SSM computation

S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','multiindex');
resonant_modes = [1 2]; % choose master spectral subspace
S.choose_E(resonant_modes);
order = 3;                  % SSM expansion order

%% numerical computation
tic
[W_0,R_0] = S.compute_whisker(order);
toc

Om = 0.98*imag(D(1));
tic
[W_1,R_1] = S.compute_perturbed_whisker(0,W_0,R_0,Om);
toc
save('multi_semi.mat','W_0','R_0','W_1','R_1');