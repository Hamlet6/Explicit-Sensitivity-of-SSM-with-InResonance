function [FRC1,FRC2,FRC3,FRC_nonauto] = ssmtool_frcs(m,b1,b2,zeta)

% build model
[mass,damp,stiff,fnl,fext] = build_model(m,zeta,b1,b2);
DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl',fnl);
set(DS.Options,'Emax',1,'Nmax',2,'notation','multiindex');
% Forcing
epsilon = 1e-2;
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(coeffs, kappas, epsilon);
% Linear Modal analysis
[V,D,W] = DS.linear_spectral_analysis();

%% SSM computation
S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','multiindex');
resonant_modes = [1 2]; % choose master spectral subspace
mFreq = 1;              % internal resonance relation vector
order = 3;                  % SSM expansion order
outdof = [1 2 3 4];           % outdof for output

set(S.FRCOptions,'sampStyle', 'cocoBD');                         % sampling style
set(S.contOptions, 'PtMX', 200, 'ItMX', 15, 'h_max', 0.1,'h_min',0.01);    % continuation setting
set(S.FRCOptions, 'nCycle',5000, 'initialSolver', 'fsolve','nt',1280);     % initial solution scheme
set(S.FRCOptions, 'coordinates', 'polar');  % two coordinate representations
% S.choose_E(resonant_modes);
% [W,R] = S.compute_whisker(order);
set(S.Options,'contribNonAuto',true);
freqRange = [1.3 1.6];
FRC1 = S.SSM_isol2ep('isol1',resonant_modes,order,mFreq,'freq',freqRange,outdof);

%%
% check convergence
FRC2 = S.SSM_isol2ep('isol2',resonant_modes,order+2,mFreq,'freq',freqRange,outdof);
FRC3 = S.SSM_isol2ep('isol3',resonant_modes,order+4,mFreq,'freq',freqRange,outdof);

ST = cell(2,1); ST{1} = {'r--','LineWidth',1.0}; ST{2} = {'r-','LineWidth',1.0};  % stable
ST2 = cell(2,1); ST2{1} = {'b--','LineWidth',1.0}; ST2{2} = {'b-','LineWidth',1.0};  % stable
ST3 = cell(2,1); ST3{1} = {'k--','LineWidth',1.0}; ST3{2} = {'k-','LineWidth',1.0};  % stable
figure; hold on;
plot_stab_lines(FRC1.om,FRC1.rho,FRC1.st,ST,'Unstable-O(3)','Stable-O(3)')
plot_stab_lines(FRC2.om,FRC2.rho,FRC2.st,ST2,'Unstable-O(5)','Stable-O(5)')
plot_stab_lines(FRC3.om,FRC3.rho,FRC3.st,ST3,'Unstable-O(7)','Stable-O(7)')
xlim(freqRange)

figure; hold on;
plot_stab_lines(FRC1.om,FRC1.Aout_frc(:,1),FRC1.st,ST,'Unstable-O(3)','Stable-O(3)')
plot_stab_lines(FRC2.om,FRC2.Aout_frc(:,1),FRC2.st,ST2,'Unstable-O(5)','Stable-O(5)')
plot_stab_lines(FRC3.om,FRC3.Aout_frc(:,1),FRC3.st,ST3,'Unstable-O(7)','Stable-O(7)')
xlim(freqRange)

% check wheather is necessary to include non-autonomous part
set(S.Options,'contribNonAuto',false);
freqRange = [1.3 1.6];
FRC_nonauto = S.SSM_isol2ep('isol-nonauto',resonant_modes,order,mFreq,'freq',freqRange,outdof);
figure; hold on;
plot_stab_lines(FRC1.om,FRC1.Aout_frc(:,1),FRC1.st,ST,'TV-Unstable','TV-Stable')
plot_stab_lines(FRC_nonauto.om,FRC_nonauto.Aout_frc(:,1),FRC_nonauto.st,ST2,'TI-Unstable','TI-Stable')
figure; hold on;
plot_stab_lines(FRC1.om,FRC1.Aout_frc(:,2),FRC1.st,ST,'TV-Unstable','TV-Stable')
plot_stab_lines(FRC_nonauto.om,FRC_nonauto.Aout_frc(:,2),FRC_nonauto.st,ST2,'TI-Unstable','TI-Stable')

end