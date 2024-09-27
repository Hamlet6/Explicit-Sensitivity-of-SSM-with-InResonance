%% Finding a 2D SSM for a 3D finite element beam
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system 
% using all phase space variables measurements. In this example, we consider the 
% clamped-clamped 3D beam.
clear all
close all


run ../../install.m
run comsol_server_ini.m
%% generate model
alpha = 0.557;
f0 = 100;
[model, M, C, K, f_0, Null, Nullf, ud, outdof, uscale] = model_mirror(alpha,f0);



% Construction of M,C,K leads to excessive data stored in Assembly class
% which in turn increases memory requirements during computation when fint
% is called -> construct fint independently

[fint,dfint] = get_fint(K,model,Null, Nullf, ud);


n = length(M);
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])
%% Dynamical system setup
% We consider the forced system
% 
% $$\mathbf{M}\ddot{\mathbf{x}}+\mathbf{C}\dot{\mathbf{x}}+\mathbf{K}\mathbf{x}+\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})=\epsilon\mathbf{f}^{ext}(\mathbf{\Omega}t),$$
% 
% which can be written in the first-order form as 
% 
% $$\mathbf{B}\dot{\mathbf{z}}	=\mathbf{Az}+\mathbf{F}(\mathbf{z})+\epsilon\mathbf{F}^{ext}(\mathbf{\phi}),\\\dot{\mathbf{\phi}}	
% =\mathbf{\Omega}$$
% 
% where
% 
% $\mathbf{z}=\left[\begin{array}{c}\mathbf{x}\\\dot{\mathbf{x}}\end{array}\right],\quad\mathbf{A}=\left[\begin{array}{cc}-\mathbf{K}  
% & \mathbf{0}\\\mathbf{0} & \mathbf{M}\end{array}\right],\mathbf{B}=\left[\begin{array}{cc}\mathbf{C}  
% & \mathbf{M}\\\mathbf{M} & \mathbf{0}\end{array}\right],\quad\quad\mathbf{F}(\mathbf{z})=\left[\begin{array}{c}\mathbf{-\mathbf{f}(\mathbf{x},\dot{\mathbf{x}})}\\\mathbf{0}\end{array}\right],\quad\mathbf{F}^{ext}(\mathbf{z},\mathbf{\phi})=\left[\begin{array}{c}\mathbf{f}^{ext}(\mathbf{\phi})\\\mathbf{0}\end{array}\right]$.

DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K);
set(DS.Options, 'Intrusion', 'none')

set(DS,'fint',fint,'dfint',dfint);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
%% 
% We assume periodic forcing of the form
% 
% $$\mathbf{f}^{ext}(\phi) = \mathbf{f}_0\cos(\phi)=\frac{\mathbf{f}_0}{2}e^{i\phi} 
% + \frac{\mathbf{f}_0}{2}e^{-i\phi}  $$
% 
% Fourier coefficients of Forcing

epsilon = 0.01;
kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
%% Linear Modal analysis and SSM setup

[V,D,W] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2];
S.choose_E(masterModes);


%%
% setup options

set(S.Options, 'reltol', 1,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',true)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 200, 'nPar', 100, 'nPsi', 200, 'rhoScale', 2 )
set(S.FRCOptions, 'method','continuation ep')
set(S.FRCOptions, 'outdof',outdof)


%% 
% choose frequency range around the first natural frequency

omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.9 1.1];
% omegaRange = [0.1, 100];
%% 
% extract forced response curve
%set(S.Options, 'PotentialNonlinearity', true)
order = [3];

FRC = S.extract_FRC('freq',omegaRange,order);
