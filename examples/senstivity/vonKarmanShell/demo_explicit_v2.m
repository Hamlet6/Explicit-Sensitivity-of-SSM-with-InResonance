%% Shallow-curved shell structure with geometric nonlinearities

% clear all; 
%% 
% *system parameters*

nDiscretization = 10; % Discretization parameter (#DOFs is proportional to the square of this number)
epsilon = 0.1;
%% generate model

[M,C,K,fnl,f_0,outdof] = build_model(nDiscretization,'type1');

fnl_han{1} = @(input) double(ttv(fnl{1},{input{1},input{2}},[2,3]));            % by tensor operation
fnl_han{2} = @(input) double(ttv(fnl{2},{input{1},input{2},input{3}},[2,3,4])); % by tensor operation

n = length(M); % number of degrees of freedom
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])
%% Dynamical system setup 
% We consider the forced system
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl_semi',fnl_han);
set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')

[V,D,W] = DS.linear_spectral_analysis();
%% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','tensor')
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
tic
[We,Re] = S.explicit_whisker();
toc

