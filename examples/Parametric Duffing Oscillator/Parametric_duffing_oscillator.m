%% Parametric Duffing Oscillator
% doi: https://doi.org/10.1006/jsvi.1995.0547
% We reproduce the setting of figure 7

clear all;clc
run ../../install.m
% change to example directory
exampleDir = fileparts(mfilename('fullpath'));
cd(exampleDir)
%% System
[alpha,beta,gamma_c,mu,q,M,C,K,Fnl,Fext] = build_model();

% For big epsilon O(1) results start to diverge at higher orders .
epsilon = 1;

%Autonomous first order mechanical system z = [x;dot(x)]
A = [-K ,0 ; 0, M];
B = [ C ,M ; M, 0];

% Dynamical System
DS = DynamicalSystem();
set(DS,'order',1)
set(DS,'A',A,'B',B,'F',Fnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')

% Analyse spectrum
[V,D,W] = DS.linear_spectral_analysis();

% External forcing
DS.add_forcing(Fext,epsilon);

%% Set up SSM object
S = SSM(DS);
set(S.Options, 'reltol', 0.5,'notation','multiindex')

%Choose Master subspace
masterModes = [1,2];
S.choose_E(masterModes);

%% Compute Autonomous SSM coefficients

order = 5; % Approximation order

%% FRC extraction routine

outdof = 1;
set(S.Options, 'reltol', 0.5,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 80, 'nPar', 20, 'nPsi', 80, 'rhoScale', 4 )
set(S.FRCOptions, 'method','level set') % 'continuation ep'
set(S.FRCOptions, 'outdof',outdof)

%% choose frequency range

omega0 = imag(S.E.spectrum(1));
OmegaRange =omega0*[0.8 1.3]; 

%% extract backbone

BB = S.extract_backbone(masterModes,OmegaRange, order);
figBB = gcf;

%% extract forced response curve

startFRCSSM = tic;
FRC = S.extract_FRC('freq',OmegaRange,order);
figFRC = gcf;
timings.FRCSSM = toc(startFRCSSM);

%% Verification: Collocation using <https://sourceforge.net/p/cocotools/wiki/Home/ coco>
% Dankowicz, H., & Schilder, F. (2013).  _Recipes for Continuation,_ SIAM Philadelphia.
% <https://doi.org/10.1137/1.9781611972573 https://doi.org/10.1137/1.9781611972573>
nCycles = 10;

coco = cocoWrapper(DS, nCycles, outdof);
set(coco,'initialGuess','forward')
set(coco.Options, 'NAdapt', 1, 'h_max', 50);
set(coco.Options,'ItMX',100,'NTST', 70,'PtMX',200); %for convergence, smaller stepsize
%and allow for more iteration subintervals and steps

figure(figFRC)
hold on
startcoco = tic;
bd3 = coco.extract_FRC(OmegaRange);
timings.cocoFRC = toc(startcoco);




