function [model, M, C, K, f_0, Null, Nullf, ud, outdof, out_full] = build_model(alpha,f0,beta)
import com.comsol.model.*
import com.comsol.model.util.*

mphopen -clear;

filename_comsol = 'beam_2d.mph';
model = mphload(filename_comsol);

model.component('comp1').physics('solid').feature('lemm1').feature('dmp1').set('alpha_dM', alpha);
model.component('comp1').physics('solid').feature('lemm1').feature('dmp1').set('beta_dK', beta);

model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);

model.study('std1').feature('stat').set('geometricNonlinearity', true);
model.study('std2').feature('time').set('geometricNonlinearity', true);


%% get MCK from time-dependent problem
MA = mphmatrix(model ,'sol2', ...
 'Out', {'K','D','E','Null','Nullf','ud','uscale'},'initmethod','sol','initsol','zero');
M_0 = MA.E; C_0 = MA.D; K_0 = MA.K; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud; uscale = MA.uscale;

M = Nullf'*M_0*Null;
% C = Nullf'*C_0*Null;
C = alpha*M;
K = Nullf'*K_0*Null;

%% set point load to get outdof

model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 1000 0]);

ML = mphmatrix(model ,'sol1','Out', {'L','K','ud','Nullf'},'initmethod','sol','initsol','zero');
L = ML.Nullf'*(ML.L-ML.K*ML.ud);
outdof = find(L);
out_full = find(ML.L);
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);
%% set boundary load to get f_ext
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 f0 0]);

MA = mphmatrix(model ,'sol1', ...
 'Out', {'K','Null','Nullf','ud','uscale','L'},'initmethod','sol','initsol','zero');
K_0 = MA.K; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud; uscale = MA.uscale;
% K = Nullf'*K_0*Null;
f_0 = Nullf'*(MA.L-MA.K*MA.ud);
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);

end