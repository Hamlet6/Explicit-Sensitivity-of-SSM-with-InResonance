%% Finding a 2D SSM for a 1D finite element beam
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system 
% using all phase space variables measurements. In this example, we consider the 
% clamped-clamped 3D beam.
clear all
close all


% run ../../install.m
% run comsol_server_ini.m

%% generate model
alpha = 0.1; % alpha = 0.5, f0 = 10; alpha = 5.85 f0 = 8500;
beta  = 1e-5;
f0 = 1000; 
[model, M, C, K, f_0, Null, Nullf, ud, outdof, out_full] = build_model(alpha,f0,beta);
% M is asymmetric for unkown reasons
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

epsilon = 0.1;
kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
% %% Linear Modal analysis and SSM setup

[V,D,W] = DS.linear_spectral_analysis();

%% factor on V and W
coefs = 1;
DS.spectrum.V = coefs*V;
DS.spectrum.W = W/coefs;

% %% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.5,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2];
S.choose_E(masterModes);

% %%
% setup options
set(S.Options, 'reltol', 1,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',false)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 300, 'nPar', 150, 'nPsi', 200, 'rhoScale', 2 )
set(S.FRCOptions, 'method','continuation ep') % 'level set' ,'continuation ep';'continuation ep', 'z0', 1e-4*[1; 1]
set(S.FRCOptions, 'outdof',outdof)


% %% 
% choose frequency range around the first natural frequency
omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.8 1.2];
set(DS.Options,'outDOF',outdof);
omegas = linspace(omegaRange(1),omegaRange(2),101);
[response,Znorm,Aout] = linear_response(DS,omegas);
figure;
plot(omegas,Aout/0.1,'k-','LineWidth',2,'DisplayName','Linear');
xlabel('\Omega'); ylabel('normalized amplitude (by thickness)')


% %% 
% extract forced response curve
%set(S.Options, 'PotentialNonlinearity', true)
order = [3 5];

FRC = S.extract_FRC('freq',omegaRange,order);
% %%
% order = [6];
% FRC2 = S.extract_FRC('freq',omegaRange,order);
%%

%% linear response
set(DS.Options,'outDOF',outdof);
omegas = [FRC{2}.Omega];
[response,Znorm,Aout] = linear_response(DS,omegas);
hold on
plot(omegas,Aout,'k-','LineWidth',2,'DisplayName','Linear')