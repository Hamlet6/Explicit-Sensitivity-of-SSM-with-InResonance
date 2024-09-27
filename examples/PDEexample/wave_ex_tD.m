%
% wave_ex.m
%
% Model exported on Oct 20 2022, 16:19 by COMSOL 5.6.0.401.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 1);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').physics.create('waeq', 'WaveEquation', 'geom1');

model.study.create('std1');
model.study('std1').create('time', 'Transient');
model.study('std1').feature('time').activate('waeq', true);

model.component('comp1').geom('geom1').run;

model.component('comp1').physics('waeq').feature('weq1').setIndex('f', 0, 0);
model.component('comp1').physics('waeq').feature('weq1').setIndex('c', '1+ux^2', 0);

model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('time').set('notlistsolnum', 1);
model.study('std1').feature('time').set('notsolnum', '1');
model.study('std1').feature('time').set('listsolnum', 1);
model.study('std1').feature('time').set('solnum', '1');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'time');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,1)');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'Default');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalvaluemethod', 'factor');
model.sol('sol1').feature('t1').set('endtimeinterpolation', true);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').attach('std1');

model.result.create('pg1', 'PlotGroup1D');
model.result('pg1').set('data', 'dset1');
model.result('pg1').create('lngr1', 'LineGraph');
model.result('pg1').feature('lngr1').set('xdata', 'expr');
model.result('pg1').feature('lngr1').set('xdataexpr', 'x');
model.result('pg1').feature('lngr1').selection.all;
model.result('pg1').feature('lngr1').set('expr', 'u');

model.sol('sol1').runAll;

model.result('pg1').run;

model.label('wave_ex.mph');

model.component('comp1').geom('geom1').create('i1', 'Interval');
model.component('comp1').geom('geom1').run('i1');
model.component('comp1').geom('geom1').run;

model.component('comp1').physics('waeq').selection.set([1]);

model.component('comp1').mesh('mesh1').automatic(false);
model.component('comp1').mesh('mesh1').autoBuildNew(false);
model.component('comp1').mesh('mesh1').feature('size').set('custom', true);

model.component('comp1').baseSystem('none');

model.component('comp1').mesh('mesh1').feature('size').set('hmax', 100);
model.component('comp1').mesh('mesh1').feature('size').set('hgrad', 1.3);
model.component('comp1').mesh('mesh1').feature('size').set('hnarrow', 1.1);
model.component('comp1').mesh('mesh1').run('size');
model.component('comp1').physics('waeq').prop('ShapeProperty').set('order', 1);
% model.component('comp1').physics('waeq').prop('ShapeProperty').set('frame', 'mesh');

model.sol('sol1').study('std1');

model.study('std1').feature('time').set('notlistsolnum', 1);
model.study('std1').feature('time').set('notsolnum', '1');
model.study('std1').feature('time').set('listsolnum', 1);
model.study('std1').feature('time').set('solnum', '1');

model.sol('sol1').feature.remove('t1');
model.sol('sol1').feature.remove('v1');
model.sol('sol1').feature.remove('st1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'time');
model.sol('sol1').create('t1', 'Time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,0.1,1)');
model.sol('sol1').feature('t1').set('plot', 'off');
model.sol('sol1').feature('t1').set('plotgroup', 'pg1');
model.sol('sol1').feature('t1').set('plotfreq', 'tout');
model.sol('sol1').feature('t1').set('probesel', 'all');
model.sol('sol1').feature('t1').set('probes', {});
model.sol('sol1').feature('t1').set('probefreq', 'tsteps');
model.sol('sol1').feature('t1').set('atolglobalvaluemethod', 'factor');
model.sol('sol1').feature('t1').set('endtimeinterpolation', true);
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('t1').feature.remove('fcDef');
model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result('pg1').run;

%%
MA = mphmatrix(model ,'sol1', ...
 'Out', {'K','L','D','E'});
M = full(MA.E); C = MA.D; K = full(MA.K);

[U,Udot]= mphgetu(model,'type','sol');
info_mesh = mphxmeshinfo(model);
info_sol = mphsolinfo(model);
[stats, meshdata] = mphmeshstats(model);
%%
% model.component('comp1').physics('waeq').feature('weq1').setIndex('c', 'ux^2', 0);
% model.physics('waeq').feature('weq1').set('ea','1');
displ  = [2;4];
u0 = [displ(1);0;0;displ(2)];
model.sol('sol1').setU(u0);
model.sol('sol1').setPVals(0);
model.sol('sol1').createSolution;

ML = mphmatrix(model, 'sol1', 'Out', {'E','K','L','Nullf','Null','ud'}, 'initmethod','sol', 'initsol', 'sol1');

Kc = ML.Nullf'*ML.K*ML.Null;
Lc = ML.Nullf'*(ML.L-ML.K*ML.ud);


fnl = -Lc - K*displ
dfnl = Kc - K
%%


