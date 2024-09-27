function [model, M, C, K, Null, Nullf, ud, outdof, out_full] = build_model(alpha)
import com.comsol.model.*
import com.comsol.model.util.*

mphopen -clear;

filename_comsol = 'arch.mph';
model = mphload(filename_comsol);

scale_unit = false;

model.component('comp1').physics('solid').feature('lemm1').feature('dmp1').set('alpha_dM', alpha);
model.component('comp1').physics('solid').feature('lemm1').feature('dmp1').set('beta_dK', 0);
%
if scale_unit
    model.component('comp1').material('mat1').propertyGroup('Enu').set('youngsmodulus', {'160e3[Pa]'});
    model.component('comp1').material('mat1').propertyGroup('def').set('density', {'2.32e-3[kg/m^3]'});
    model.component('comp1').geom('geom1').lengthUnit('m');
else
    model.component('comp1').material('mat1').propertyGroup('Enu').set('youngsmodulus', {'160e9[Pa]'});
    model.component('comp1').material('mat1').propertyGroup('def').set('density', {'2320[kg/m^3]'});
end

model.study('std1').feature('stat').set('geometricNonlinearity', true);
model.study('std2').feature('time').set('geometricNonlinearity', true);


%% get MCK from time-dependent problem
MA = mphmatrix(model ,'sol2', ...
 'Out', {'K','D','E','Null','Nullf','ud','uscale'},'initmethod','sol','initsol','zero');
M_0 = MA.E; C_0 = MA.D; K_0 = MA.K; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud; uscale = MA.uscale;

M = Nullf'*M_0*Null;
C = Nullf'*C_0*Null;
% C = alpha*M;
K = Nullf'*K_0*Null;

%% set point load to get outdof

model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 1000 0]);
% 
ML = mphmatrix(model ,'sol1','Out', {'L','K','ud','Nullf'},'initmethod','sol','initsol','zero');
L = ML.Nullf'*(ML.L-ML.K*ML.ud);
outdof = find(L);
out_full = find(ML.L);
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);

end