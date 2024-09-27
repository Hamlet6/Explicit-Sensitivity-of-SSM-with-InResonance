function [model, M, C, K, f_0, Null, Nullf, ud, outdof, out_full] = build_model(alpha,f0,beta)
import com.comsol.model.*
import com.comsol.model.util.*

mphopen -clear;

filename_comsol = 'model_3d.mph';
model = mphload(filename_comsol);

% model.component('comp1').physics('solid').feature('lemm1').feature('dmp1').set('alpha_dM', alpha);
% model.component('comp1').physics('solid').feature('lemm1').feature('dmp1').set('beta_dK', beta);

model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);

model.study('std1').feature('stat').set('geometricNonlinearity', true);
% model.study('std2').feature('time').set('geometricNonlinearity', true);


%% get MCK from time-dependent problem
MA = mphmatrix(model ,'sol2', ...
 'Out', {'K','D','E','Null','Nullf','ud','uscale'},'initmethod','sol','initsol','zero');
M_0 = MA.E; C_0 = MA.D; K_0 = MA.K; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud; uscale_dyn = MA.uscale;

M = Nullf'*M_0*Null;
% C = Nullf'*C_0*Null;
K = Nullf'*K_0*Null;
C = alpha*M+beta*K;

%% set point load to get outdof

model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 1000]);

ML = mphmatrix(model ,'sol1','Out', {'L','K','ud','Nullf'},'initmethod','sol','initsol','zero');
L = ML.Nullf'*(ML.L-ML.K*ML.ud);
outdof = find(L);
out_full = find(ML.L);

%% set boundary load to get f_ext
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 f0]);

MA = mphmatrix(model ,'sol1', ...
 'Out', {'K','Null','Nullf','ud','uscale','L'},'initmethod','sol','initsol','zero');
K_0 = MA.K; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud; uscale_sta = MA.uscale;
Ks = Nullf'*K_0*Null;
% we check whether uscale_dyn is consistent with uscale_sta and also the
% consistency of stiffness matrix
res1 = norm(uscale_dyn-uscale_sta);
fprintf('residual for scaling factor is %d\n', res1);
res2 = K-Ks;
fprintf('residual for stiffness matrix is %d\n', norm(res2(:),'inf'));
f_0 = Nullf'*(MA.L-MA.K*MA.ud);
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);

end