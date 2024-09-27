%% Arch example
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system 
% using all phase space variables measurements. In this example, we consider the 
% shallow arch geometry. Mesh examples are also included.
clear all
close all

run ../../install.m
run comsol_server_ini.m

%% generate model

omega_1 = 0.92227;
alpha = omega_1/500;

[model, M, C, K, Null, Nullf, ud, outdof, out_full] = build_model(alpha);

[fint,dfint] = get_fint(K,model,Null, Nullf, ud);

n = length(M);
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])
% %% Dynamical system setup
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K);
set(DS.Options, 'Intrusion', 'none')

set(DS,'fint',fint,'dfint',dfint);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
%%
[V,D,W] = DS.linear_spectral_analysis();
%%
masterModes = [9 10];
Phi_1 = V(1:n,masterModes(1));
f_0 = M*Phi_1;
[~, outdof] = max(Phi_1);
%%
% epsilon = 0.06*(0.92097e-0)^2;
epsilon = 0.09*omega_1^2;

kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
% %% Linear Modal analysis and SSM setup


% %%
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')

S.choose_E(masterModes);


% %%
% setup options

set(S.Options, 'reltol', 1,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',false)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 1000, 'nPar', 500, 'nPsi', 200, 'rhoScale', 2 )
set(S.FRCOptions, 'method', 'level set') % 'level set' ,'continuation ep';'continuation ep', 'z0', 1e-4*[1; 1]
set(S.FRCOptions, 'outdof',outdof)
set(DS.Options, 'outDOF',outdof)

% %% 
% choose frequency range around the first natural frequency

omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.95 1.05];
% %% 
% extract forced response curve
%set(S.Options, 'PotentialNonlinearity', true)
order = [3,5];
%% Create linear response analysis
omegas = linspace(omegaRange(1),omegaRange(2), 500);
[response,Znorm,Aout,Aa] = linear_response(DS,omegas,Phi_1);
%%
FRC = S.extract_FRC('freq',omegaRange,order);

%% plot FRC
colors = get(0,'defaultaxescolororder');
f1 = figure('Name','Norm');
f2 = figure('Name',['Amplitude at DOFs ' num2str(outdof)]);
figs = [f1, f2];

plot_FRC(FRC{1,1},outdof,3,'freq','circles', figs, colors(1,:))
%% plot linear
figure; hold on
stab_plot(omegas/omega0,Aout/6.4,ones(500,1),3,colors(1,:));
