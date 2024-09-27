%% MYMEMS
clearvars
% We compute SSM compute for a geometrically nonlinear Timoshenko Beam
order = [3];
clc
% Generate model
[M,C,K,fnltens,fext,outdof] = build_model('z');
n = length(M);


mode = 1;
masterModes = [2*mode-1,2*mode]; 
epsilon = 5e-6; % 

% Dynamical system setup
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnltens);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')

% Linear Modal analysis and SSM setup

[V,D,W] = DS.linear_spectral_analysis();

kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);

%% 
% *Choose Master subspace*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
S.choose_E(masterModes);
%% 
% Setup options

set(S.Options, 'reltol', 0.5,'IRtol',0.05,'notation', 'multiindex','contribNonAuto',false)
set(S.FRCOptions, 'nt', 2^8, 'nRho', 200, 'nPar', 200, 'nPsi', 100, 'rhoScale',1)
set(S.FRCOptions, 'outdof',outdof)
set(S.FRCOptions, 'method','level set') 

%% 
% Choose frequency range

omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.9,1.1];
%% 
% extract forced response curve

FRC = S.extract_FRC('freq',omegaRange,order);


