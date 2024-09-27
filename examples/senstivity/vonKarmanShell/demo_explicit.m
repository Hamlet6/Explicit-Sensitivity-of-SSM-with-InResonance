%% Shallow-curved shell structure with geometric nonlinearities

% clear all; 
%% 
% *system parameters*

nDiscretization = 10; % Discretization parameter (#DOFs is proportional to the square of this number)
epsilon = 0.1;
%% generate model

[M,C,K,fnl,dfnl,f_0,outdof] = build_model_functional(nDiscretization,'type1');
n = length(M); % number of degrees of freedom
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])
%% Dynamical system setup 
% We consider the forced system
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl_semi',fnl,'dfnl_semi',dfnl);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
set(DS.Options,'Intrusion','semi')

[V,D,W] = DS.linear_spectral_analysis();

% recale Vs and Ws such that phi'Mphi = I
phi = DS.spectrum.V(1:n,1);
tmp = sqrt(phi'*M*phi);
DS.spectrum.V = V/tmp;
DS.spectrum.W = W*tmp;

%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex','solver','backslash')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [1,2]; 
S.choose_E(masterModes);
%% compute SSMs
% Obtaining *forced response curve* in reduced-polar coordinate

order = 3; % Approximation order
%% numerical computation
tic
[W_0,R_0] = S.compute_whisker(order);
toc
% 
% save('multi-semi.mat','W_0','R_0')

%% explicit computation
[We,Re] = S.explicit_whisker();
