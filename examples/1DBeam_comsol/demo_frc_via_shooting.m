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

set(DS,'fnl_non',fint,'dfnl_non',dfint);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')

epsilon = 0.1;
kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
% %% Linear Modal analysis and SSM setup

[V,D,W] = DS.linear_spectral_analysis();

%% linear FRC
% choose frequency range around the first natural frequency
omega0 = imag(D(1));
omegaRange = omega0*[0.8 1.2];
set(DS.Options,'outDOF',outdof);
omegas = linspace(omegaRange(1),omegaRange(2),101);
[response,Znorm,Aout] = linear_response(DS,omegas);
figure;
plot(omegas,Aout/0.1,'k-','LineWidth',2,'DisplayName','Linear');
xlabel('\Omega'); ylabel('normalized amplitude (by thickness)')

plot(omegas,Aout,'k-','LineWidth',2,'DisplayName','Linear')

%% nonlinear FRC via shooting
nCycles = 10;
coco = cocoWrapper(DS, nCycles, outdof);
set(coco.Options, 'PtMX', 100);
set(coco.Options, 'h_max', 20, 'NSV',10);
set(coco, 'initialGuess', 'linear');
set(coco.Options, 'Nsteps', 50, 'RelTol', 1e-5);
set(coco.Options, 'dir_name', 'ep01');
startcoco = tic;
coco.forward_FRC(omegaRange);
timings.cocoFRC = toc(startcoco);

hold on;
coco_plot_bd('ep01.forwardFRC','omega','amp5');
