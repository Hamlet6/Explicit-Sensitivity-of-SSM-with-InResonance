function [alpha,beta,gamma_c,mu,q,M,C,K,Fnl,Fext] = build_model()
%% Parameter Versions from the paper
%% With forcing figure 7
[alpha,beta,gamma_c,mu,q] = deal(0.1,0.05,0.1,0.2,0.2);

%% Dynamical System

M = 1;
C = -alpha;
K = 1;

%% Nonlinearity
Fnl(3).coeffs = -[beta,gamma_c;0,0];
Fnl(3).ind    = [0,3;3,0];


%% External forcing
Fext.data = set_forcing_coefficients(gamma_c,mu,q);

end


function [data] = set_forcing_coefficients(gamma,mu,q)
% Set external forcing coefficients
% Multi-indices are stored in rows, coefficients in columns

f_ord = 4; %nonlinearity x^3 including zeroth order


idle = repmat(struct('coeffs',[],'ind',[]),f_ord,1);
data = repmat(struct('kappa',[],'F_n_k',idle),4,1);

% kappa_1
data(1).kappa = 1;
data(1).F_n_k(1).coeffs = [q/2;0];
data(1).F_n_k(1).ind    = [0,0];

% kappa_2
data(2).kappa = -1;
data(2).F_n_k(1).coeffs = [q/2;0];
data(2).F_n_k(1).ind    = [0,0];

% kappa_3
data(3).kappa = 2;
data(3).F_n_k(2).coeffs = [mu/2;0];
data(3).F_n_k(2).ind    = [1,0];

data(3).F_n_k(4).coeffs = [gamma*mu/2;0];
data(3).F_n_k(4).ind    = [3,0];

% kappa_4
data(4).kappa = -2;
data(4).F_n_k(2).coeffs = [mu/2;0];
data(4).F_n_k(2).ind    = [1,0];

data(4).F_n_k(4).coeffs = [gamma*mu/2;0];
data(4).F_n_k(4).ind    = [3,0];
end
