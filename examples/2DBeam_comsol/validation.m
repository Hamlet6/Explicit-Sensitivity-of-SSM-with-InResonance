%% Finding a 2D SSM for a 2D finite element beam
% This is an example of how to reconstruct a slow 2D SSM of a mechanical system 
% using all phase space variables measurements. In this example, we consider the 
% clamped-clamped 3D beam.
clear all
close all


run ../../install.m
run comsol_server_ini.m

%% generate model
alpha = 5.85; % alpha = 0.5, f0 = 10; alpha = 5.85 f0 = 8500;
beta = 0;
f0 = 7500; % 600 11000
[model, M, C, K, f_0, Null, Nullf, ud, outdof, out_full] = build_model(alpha,f0,beta);

[fint,dfint] = get_fint(K,model,Null, Nullf, ud);

n = length(M);
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])
% %% Dynamical system setup
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K);
set(DS.Options, 'Intrusion', 'none')

set(DS,'fint',fint,'dfint',dfint);
set(DS.Options,'Emax',5,'Nmax',10,'notation','multiindex')
% set(DS.Options,'Emax',5,'Nmax',10,'notation','tensor')
% %%
% f_0(outdof) = 300; % 28000 % alpha = 0.25- f_0 = 300;

epsilon = 0.01;
kappas = [-1; 1];
coeffs = [f_0 f_0]/2;
DS.add_forcing(coeffs, kappas,epsilon);
% %% Linear Modal analysis and SSM setup

[V,D,W] = DS.linear_spectral_analysis();
% %% 
% *Choose Master subspace (perform resonance analysis)*

S = SSM(DS);
set(S.Options, 'reltol', 0.1,'notation','multiindex')
% set(S.Options, 'reltol', 0.1,'notation','tensor')
masterModes = [7,8]; % 7 8
S.choose_E(masterModes);


% %%
% setup options

set(S.Options, 'reltol', 1,'IRtol',0.02,'notation', 'multiindex','contribNonAuto',false)
set(S.FRCOptions, 'nt', 2^7, 'nRho', 300, 'nPar', 150, 'nPsi', 200, 'rhoScale', 2 )
set(S.FRCOptions, 'method','level set') % 'level set' ,'continuation ep';'continuation ep', 'z0', 1e-4*[1; 1]
set(S.FRCOptions, 'outdof',outdof)


% %% 
% choose frequency range around the first natural frequency

omega0 = imag(S.E.spectrum(1));
omegaRange = omega0*[0.8 1.2];
% %% 
% extract forced response curve
%set(S.Options, 'PotentialNonlinearity', true)
order = [3,5];

FRC = S.extract_FRC('freq',omegaRange,order);
% %%
% order = [6];
% FRC2 = S.extract_FRC('freq',omegaRange,order);
%%
colors = get(0,'defaultaxescolororder');
f1 = figure('Name','Norm');
f2 = figure('Name',['Amplitude at DOFs ' num2str(outdof)]);
figs = [f1, f2];

plot_FRC(FRC{1,1},outdof,order,'freq','circles', figs, colors(1,:))
hold on
plot_FRC(FRC{1,2},outdof,order,'freq','circles', figs, colors(2,:))
% hold on
% plot_FRC(FRC{1,1},outdof,order,'freq','lines', figs, colors(3,:))
% hold on
% plot_FRC(FRC{1,2},outdof,order,'freq','lines', figs, colors(4,:))
%%
idx_set = [72 74];
u_outi = [];
for i = 1:length(idx_set)
    if i == 1
        created = 0;
    else
        created = 1;
    end
idx = idx_set(i);
% created =  0;
if created == 0
 model.component('comp1').func.create('an1', 'Analytic');
end
ic = FRC{1, 1}(idx).Zic;
freq = FRC{1, 1}(idx).Omega;
model.component('comp1').func('an1').set('expr', ['cos(',num2str(freq),'*x)']);
amp = num2str(f0*epsilon);
model.component('comp1').physics('solid').feature('pl1').set('Fp', {'0', [amp,'*an1(t)'], '0'});
nspc = 50;  % number of samples per (forcing) cycle
nc = 300;    % number of cycles
rtol = 1e-6; % relative tolerance for time integration

%% Parameters

% time simulation parameters
f1 = freq/2/pi;
model.param.set('f1', [num2str(f1) ' [Hz]'], 'forcing frequency');
model.param.set('Tc', '1/f1',       'cycle time');
model.param.set('nspc', nspc,       'number of samples per cycle');
model.param.set('dt', 'Tc/nspc',    'time step');                   
model.param.set('nc', nc,           'number of cycles');
model.param.set('Tend', 'Tc*nc',    'transient end time');
Tc = 1/f1;
dt = Tc/nspc;
Tend = Tc*nc;
%%
if created == 0
    model.study.create('std3');
    model.study('std3').create('time', 'Transient');
end
model.study('std3').feature('time').activate('solid', true);

model.study("std3").run(); % time transient with 1 time step (used only to define sol1)
u0 = Null*ic(1:n)'+ud;
ud0 = Null*ic(n+1:end)'+ud;
% modify sol1 using our ICs
U = u0;
UDOT = ud0;
% U = mphgetu(model,'soltag','sol1');
model.sol('sol3').setU(U);
model.sol('sol3').setUDot(UDOT);
model.sol('sol3').setPNames('t');
model.sol('sol3').setPVals(1);
model.sol('sol3').createSolution;
%%
NL = 'on';

if created == 0
    model.study.create('std4');
    model.study('std4').create('time', 'Transient');
end
model.study('std4').feature('time').set('tlist', 'range(0,0.01,2)');
model.study('std4').feature('time').set('geometricNonlinearity', NL);
model.study('std4').feature('time').set('probefreq', 'tout');

model.study('std4').feature('time').set('useinitsol', true);   % use initial solution
model.study('std4').feature('time').set('initmethod', 'sol');  % by selecting a previous solution
model.study('std4').feature('time').set('initstudy', 'std3');  % which is std1 (study 1)
model.study('std4').feature('time').set('solnum', 'first');    % selecting the first solution

% SOLUTION 4
if created == 0
model.sol.create('sol4');
model.sol('sol4').study('std4');
model.sol('sol4').attach('std4');
model.sol('sol4').create('st1', 'StudyStep');
model.sol('sol4').create('v1', 'Variables');
model.sol('sol4').create('t1', 'Time');
model.sol('sol4').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol4').feature('t1').create('d1', 'Direct');
model.sol('sol4').feature('t1').create('i1', 'Iterative');
model.sol('sol4').feature('t1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol4').feature('t1').feature('i1').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol4').feature('t1').feature('i1').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol4').feature('t1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol4').feature('t1').feature.remove('fcDef');
end
model.sol('sol4').attach('std4');
model.sol('sol4').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol4').feature('v1').label('Dependent Variables 1.1');
model.sol('sol4').feature('t1').set('probefreq', 'tout');
model.sol("sol4").feature("t1").set("timemethod", "genalpha"); % solver
model.sol("sol4").feature("t1").set("tstepsgenalpha", "strict"); % force the solver to use the time steps defined by the user

model.sol("sol4").feature("t1").set("initialstepgenalphaactive", true);
model.sol("sol4").feature("t1").set("initialstepgenalpha", "dt");

model.study("std4").feature("time").set("usertol", true);
model.study("std4").feature("time").set("rtol", rtol);
%% 
model.study('std4').feature('time').set('tlist', 'range(0,dt,Tend)');
%%
tic
model.study("std4").run();
toc
% %%

sol4info = mphsolinfo(model,'soltag','sol4','NU','on');
NU = sol4info.sizesolvals;
[u, udot] = mphgetu(model,'soltag','sol4','solnum',1:1:NU);
u_out = u(out_full,:);
udot_out = udot(out_full,:);
u_outi = [u_outi;u_out;udot_out];
% %%
figure
plot(0:dt:Tend, u_out, 'linewidth', 1)
xlabel('t  [sec]')
ylabel('u [m]')
title(['time series, IC'])
axis tight

plotXnt(u_out,nspc)
disp(['i = ',num2str(i),';'])
% amp = [amp,max(u_out(end-1000:end))];
end
ampl = [];
omegas = [];
for j = 1:length(idx_set)
    ampl = [ampl,max(u_outi(j*2-1,end-1000:end))];
    omegas = [omegas,FRC{1,1}(idx_set(j)).Omega];
end
%%
% save('al5_85f07500/idx.mat','u_outi','alpha','f0','idx_set','FRC')
%%
% figure
% plot(0:dt:Tend, u_outi(5,:), 'linewidth', 1)
% xlabel('t  [sec]')
% ylabel('u [m]')
% title(['time series, IC'])
% axis tight
%% plot together with full simulation results and linear response
colors = get(0,'defaultaxescolororder');
f1 = figure('Name','Norm');
f2 = figure('Name',['Amplitude at DOFs ' num2str(outdof)]);
figs = [f1, f2];
plot_FRC(FRC{1,1},outdof,3,'freq','circles', figs, colors(1,:))
hold on
plot_FRC(FRC{1,2},outdof,5,'freq','circles', figs, colors(2,:))
hold on
plot_FRC(FRCLinear{1,1},outdof,3,'freq','lines', figs, colors(3,:))
hold on
plot(omegas,ampl,'sqr','color','black','DisplayName','full simulation')
%% validation of nonlinear stiffness & full system simulation (time dependent)
f0 = 7500;

NL = 'on';
rtol = 1e-6;
model.param.set('dt', '0.1',    'time step');                   
model.param.set('Tend', '3',    'transient end time');

model.component('comp1').physics('solid').feature('pl1').set('Fp', {'0', [num2str(f0*epsilon)], '0'});
created = 0;

if created == 0
model.study.create('std5');
model.study('std5').create('time', 'Transient');
end
model.study('std5').feature('time').set('tlist', 'range(0,dt,Tend)');
model.study('std5').feature('time').set('geometricNonlinearity', NL);
model.study('std5').feature('time').set('probefreq', 'tout');

% model.study('std5').feature('time').set('useinitsol', true);   % use initial solution
% model.study('std5').feature('time').set('initmethod', 'sol');  % by selecting a previous solution
% model.study('std5').feature('time').set('initstudy', 'std3');  % which is std1 (study 1)
% model.study('std5').feature('time').set('solnum', 'first');    % selecting the first solution

% SOLUTION 5
if created == 0
model.sol.create('sol5');
model.sol('sol5').study('std5');
model.sol('sol5').attach('std5');
model.sol('sol5').create('st1', 'StudyStep');
model.sol('sol5').create('v1', 'Variables');
model.sol('sol5').create('t1', 'Time');
model.sol('sol5').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol5').feature('t1').create('d1', 'Direct');
model.sol('sol5').feature('t1').create('i1', 'Iterative');
model.sol('sol5').feature('t1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol5').feature('t1').feature('i1').feature('mg1').feature('pr').create('so1', 'SOR');
model.sol('sol5').feature('t1').feature('i1').feature('mg1').feature('po').create('so1', 'SOR');
model.sol('sol5').feature('t1').feature('i1').feature('mg1').feature('cs').create('d1', 'Direct');
model.sol('sol5').feature('t1').feature.remove('fcDef');
end
model.sol('sol5').attach('std5');
model.sol('sol5').feature('st1').label('Compile Equations: Time Dependent');
model.sol('sol5').feature('v1').label('Dependent Variables 1.1');
model.sol('sol5').feature('t1').set('probefreq', 'tout');
model.sol("sol5").feature("t1").set("timemethod", "genalpha"); % solver
model.sol("sol5").feature("t1").set("tstepsgenalpha", "strict"); % force the solver to use the time steps defined by the user

model.sol("sol5").feature("t1").set("initialstepgenalphaactive", true);
model.sol("sol5").feature("t1").set("initialstepgenalpha", "dt");

model.study("std5").feature("time").set("usertol", true);
model.study("std5").feature("time").set("rtol", rtol);

tic
model.study("std5").run();
toc
% %%

sol5info = mphsolinfo(model,'soltag','sol5','NU','on');
NU5 = sol5info.sizesolvals;
[u_sim, udot_sim] = mphgetu(model,'soltag','sol5','solnum',1:NU5);

% %%
figure
plot(0:0.1:3,u_sim(out_full,:), 'linewidth', 1)
xlabel('t  [sec]')
ylabel('u [m]')
title(['time series, IC'])
axis tight
%% calculate by stationary solver
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 f0*epsilon 0]);

created = 0;
if created == 0
model.study.create('std6');
model.study('std6').create('stat', 'Stationary');
model.study('std6').feature('stat').activate('solid', true);
model.sol.create('sol6');
end
model.sol('sol6').study('std6');

model.study('std6').feature('stat').set('notlistsolnum', 1);
model.study('std6').feature('stat').set('notsolnum', '1');
model.study('std6').feature('stat').set('listsolnum', 1);
model.study('std6').feature('stat').set('solnum', '1');
if created == 0
model.sol('sol6').create('st1', 'StudyStep');
model.sol('sol6').feature('st1').set('study', 'std6');
model.sol('sol6').feature('st1').set('studystep', 'stat');
model.sol('sol6').create('v1', 'Variables');
model.sol('sol6').feature('v1').set('control', 'stat');
model.sol('sol6').create('s1', 'Stationary');
model.sol('sol6').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol6').feature('s1').create('fc1', 'FullyCoupled');
end
model.study('std6').feature('stat').set('geometricNonlinearity', true);

model.sol('sol6').feature('s1').feature('fc1').set('termonres', 'auto');
model.sol('sol6').feature('s1').feature('fc1').set('reserrfact', 1000);
model.sol('sol6').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol6').feature('s1').feature('fc1').set('termonres', 'auto');
model.sol('sol6').feature('s1').feature('fc1').set('reserrfact', 1000);
model.sol('sol6').feature('s1').feature.remove('fcDef');
model.sol('sol6').attach('std6');

tic
model.study("std6").run();
toc
% %%
% model.sol('sol6').runAll;

[u_stationary] = mphgetu(model,'soltag','sol6');

%% calculate by Newton method
model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);
x0 = zeros(n,1);
f0 = 7500;
f_0(outdof) = f0*epsilon;
func = @(x) K*x - f_0;
func_nonlinear =  @(x) K*x + fint(x) - f_0;
% u + ud;
u_linear = Null*newton_raph(func, x0, dfint, K) + ud;
u_Nlinear = Null*newton_raph(func_nonlinear, x0, dfint, K) + ud;
% u_sim(:,end);
%% relative error
diff_norm = norm(u_linear - u_Nlinear)
% diff_norm1 = norm(u_sim(:,end) - u_linear)
%%
figure
plot(u_N,'DisplayName','Newton')
hold on
plot(u_sim(:,end),'DisplayName','comsol')
hold on
plot(u_stationary(:,end),'DisplayName','comsol stationary')
legend

