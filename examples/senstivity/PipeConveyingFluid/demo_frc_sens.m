clear all;

nmodes = 4;
n = nmodes;
phiend = zeros(n,1);
lamda  = zeros(n,1);
lamda(1)=1.8751040687119611664453082410782141625701117335311;
lamda(2)=4.6940911329741745764363917780198120493898967375458;
lamda(3)=7.8547574382376125648610085827645704578485419292300;
lamda(4)=10.995540734875466990667349107854702939612972774652;
lamda(5)=14.137168391046470580917046812551772068603076792975;
lamda(6)=17.278759532088236333543928414375822085934519635550;
lamda(7)=20.420352251041250994415811947947837046137288894544;
lamda(8)=23.561944901806443501520253240198075517031265990051;
for k=1:n
   x    = 1;
   phin = cos(lamda(k)*x)-cosh(lamda(k)*x)-(cos(lamda(k))+cosh(lamda(k)))/(sin(lamda(k))+sinh(lamda(k)))*...
        (sin(lamda(k)*x)-sinh(lamda(k)*x));%clamped-free
   phiend(k) = phin;
end
Q = phiend*phiend';
% disp at the free end
obsfun = @(x,mapx,nmodes) mapx*x(1:nmodes,:);
obs = @(x) obsfun(x,phiend',n);
%% 
% Create dynamical model
nmodes = 4;
fcload = 1;
flowspeed = 6;
miu = 0;
Gamma=0;
alpha=0.001;
beta=0.2;
[M, C, K, fnl, fext] = build_model(nmodes,flowspeed,beta,miu,Gamma,alpha,'clamped-free','nonlinear_damp');
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl',fnl);
set(DS.Options,'Emax',10,'Nmax',10,'notation','multiindex')
set(DS.Options,'RayleighDamping',false)
[V,D,W] = DS.linear_spectral_analysis();

Omega = 13;
epsilon = 0.003;
kappas = [-1; 1];
coeffs = [fext fext]/2;
DS.add_forcing(Omega^2*coeffs, kappas, epsilon);

% SSM computation
S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','multiindex');
set(S.Options, 'contribNonAuto',true);
resonant_modes = [1 2]; % choose master spectral subspace
mFreq = 1;              % internal resonance relation vector
order = 3;              % SSM expansion order
set(S.FRCOptions,'sampStyle', 'cocoBD', 'nt', 1280);                       % sampling style
set(S.contOptions, 'PtMX', 200, 'ItMX', 15, 'h_max', 0.1,'h_min',0.01);    % continuation setting
set(S.FRCOptions, 'nCycle',5000, 'initialSolver', 'fsolve');     % initial solution scheme
set(S.FRCOptions, 'coordinates', 'polar');  % two coordinate representations
% S.choose_E(resonant_modes);
% [W,R] = S.compute_whisker(order);
freqRange = [13 15]; outdof = 1:4;
FRC = S.SSM_isol2ep('isol',resonant_modes,order,mFreq,'freq',freqRange,outdof);

%% explicit computation
[M, C, K, fnl,fext,DM,DC,DK,dfnl2,dfnl3,dfext] = build_model_semi(nmodes,flowspeed,beta,alpha,'clamped-free');
DS = DynamicalSystem();
set(DS,'M',M,'C',C,'K',K,'fnl_semi',fnl,'nl_damp',true);
set(DS.Options,'Emax',10,'Nmax',10,'notation','multiindex');
set(DS.Options,'RayleighDamping',false);
set(DS,'fext',Omega^2*fext,'Omega',Omega,'dfext',dfext);
S = SSM(DS);
set(S.Options, 'reltol', 1,'notation','multiindex');
set(S.Options, 'contribNonAuto',true);
resonant_modes = [1 2]; % choose master spectral subspace
S.choose_E(resonant_modes);
tic
[We,Re] = S.explicit_whisker();
toc

% SSM senstivity computation
set(DS,'DM',DM,'DC',DC,'DK',DK);
set(DS,'Dfnl2_semi',dfnl2,'Dfnl3_semi',dfnl3);
% set(S.Options,'solver','backslash');
[DW,DR] = S.explicit_senstivity_whisker(We,Re);

optdof = 1:nmodes; 
% optdof = outdof;
[frce,dfrce] = S.po_sensitivity(We,Re,DW,DR,optdof,Q,epsilon);

%% validation via central finite difference
%% w.r.t. u
uE  = DS.spectrum.W(:,1);
psi = conj(uE(1:DS.n)); lo = (2*DS.spectrum.Lambda(1)*M+C).'*psi;
%%
du  = 0.001*flowspeed; % 0.001 is default
frc_up = comp_frc(nmodes,flowspeed+du,beta,alpha,lo,Omega,Q,epsilon);
frc_um = comp_frc(nmodes,flowspeed-du,beta,alpha,lo,Omega,Q,epsilon);
drho_du = (frc_up.rho - frc_um.rho)/(2*du);
dth_du  = (frc_up.th - frc_um.th)/(2*du);
damp_du = (frc_up.amp - frc_um.amp)/(2*du);
drhomax_du = (frc_up.rhomax-frc_um.rhomax)/(2*du);
dommax_du  = (frc_up.Ommax-frc_um.Ommax)/(2*du);
% calculate derivatives
fprintf('relative error for drho_du is %d\n', abs((dfrce.Drho(1)-drho_du)/dfrce.Drho(1))); 
fprintf('relative error for dth_du is %d\n', abs((dfrce.Dth(1)-dth_du)/dfrce.Dth(1))); 
fprintf('relative error for damp_du is %d\n', abs((dfrce.Damp(1)-damp_du)/dfrce.Damp(1)));
fprintf('realtive error for drhomax_du is %d\n', abs((dfrce.Drhomax(1)-drhomax_du)/dfrce.Drhomax(1)));
fprintf('relative error for dommax_du is %d\n', abs((dfrce.DOmmax(1)-dommax_du)/dfrce.DOmmax(1)));
% NOTE!!! The last one converges very slowly
%% w.r.t beta
dbe = 0.001*beta; 
frc_bep = comp_frc(nmodes,flowspeed,beta+dbe,alpha,lo,Omega,Q,epsilon);
frc_bem = comp_frc(nmodes,flowspeed,beta-dbe,alpha,lo,Omega,Q,epsilon);
drho_dbe = (frc_bep.rho - frc_bem.rho)/(2*dbe);
dth_dbe  = (frc_bep.th - frc_bem.th)/(2*dbe);
damp_dbe = (frc_bep.amp - frc_bem.amp)/(2*dbe);
drhomax_dbe = (frc_bep.rhomax-frc_bem.rhomax)/(2*dbe);
dommax_dbe  = (frc_bep.Ommax-frc_bem.Ommax)/(2*dbe);
% calculate derivatives
fprintf('relative error for drho_dbe is %d\n', abs((dfrce.Drho(2)-drho_dbe)/dfrce.Drho(2))); 
fprintf('relative error for dth_dbe is %d\n', abs((dfrce.Dth(2)-dth_dbe)/dfrce.Dth(2))); 
fprintf('relative error for damp_dbe is %d\n', abs((dfrce.Damp(2)-damp_dbe)/dfrce.Damp(2)));
fprintf('realtive error for drhomax_dbe is %d\n', abs((dfrce.Drhomax(2)-drhomax_dbe)/dfrce.Drhomax(2)));
fprintf('relative error for dommax_dbe is %d\n', abs((dfrce.DOmmax(2)-dommax_dbe)/dfrce.DOmmax(2)));

%% w.r.t alpha
dal = 0.01*alpha; 
frc_alp = comp_frc(nmodes,flowspeed,beta,alpha+dal,lo,Omega,Q,epsilon);
frc_alm = comp_frc(nmodes,flowspeed,beta,alpha-dal,lo,Omega,Q,epsilon);
drho_dal = (frc_alp.rho - frc_alm.rho)/(2*dal);
dth_dal  = (frc_alp.th - frc_alm.th)/(2*dal);
damp_dal = (frc_alp.amp - frc_alm.amp)/(2*dal);
drhomax_dal = (frc_alp.rhomax-frc_alm.rhomax)/(2*dal);
dommax_dal  = (frc_alp.Ommax-frc_alm.Ommax)/(2*dal);
% calculate derivatives
fprintf('relative error for drho_dal is %d\n', abs((dfrce.Drho(3)-drho_dal)/dfrce.Drho(3))); 
fprintf('relative error for dth_dal is %d\n', abs((dfrce.Dth(3)-dth_dal)/dfrce.Dth(3))); 
fprintf('relative error for damp_dal is %d\n', abs((dfrce.Damp(3)-damp_dal)/dfrce.Damp(3)));
fprintf('realtive error for drhomax_dal is %d\n', abs((dfrce.Drhomax(3)-drhomax_dal)/dfrce.Drhomax(3)));
fprintf('relative error for dommax_dal is %d\n', abs((dfrce.DOmmax(3)-dommax_dal)/dfrce.DOmmax(3)));

%% perturbed FRC
pidx = 2;
pfracs = [-0.002 0 0.002];
pFRCs = S.FRC_sensitivity(We,Re,DW,DR,optdof,Q,epsilon,freqRange,pidx,pfracs,'isola');

ST1 = cell(2,1);
ST1{1} = {'r--','LineWidth',1.0}; % unstable
ST1{2} = {'r-','LineWidth',1.0};  % stable
ST2 = cell(2,1);
ST2{1} = {'b--','LineWidth',1.0}; % unstable
ST2{2} = {'b-','LineWidth',1.0};  % stable
ST3 = cell(2,1);
ST3{1} = {'k--','LineWidth',1.0}; % unstable
ST3{2} = {'k-','LineWidth',1.0};  % stable

figure; hold on
plot_stab_lines(pFRCs{1}.om,pFRCs{1}.rho,pFRCs{1}.st,ST1);
plot_stab_lines(pFRCs{2}.om,pFRCs{2}.rho,pFRCs{2}.st,ST2);
plot_stab_lines(pFRCs{3}.om,pFRCs{3}.rho,pFRCs{3}.st,ST3);
xlim(freqRange); grid on; box on
xlabel('\Omega'); ylabel('\rho');

figure; hold on
plot_stab_lines(pFRCs{1}.om,pFRCs{1}.ampL2,pFRCs{1}.st,ST1);
plot_stab_lines(pFRCs{2}.om,pFRCs{2}.ampL2,pFRCs{2}.st,ST2);
plot_stab_lines(pFRCs{3}.om,pFRCs{3}.ampL2,pFRCs{3}.st,ST3);
xlim(freqRange); grid on; box on
xlabel('\Omega'); ylabel('Amp');

% direct computation
[~,frcp1] = comp_frc(n,flowspeed,beta+pfracs(1),alpha,lo,Omega,Q,epsilon,freqRange,'isola');
[~,frcp3] = comp_frc(n,flowspeed,beta+pfracs(3),alpha,lo,Omega,Q,epsilon,freqRange,'isola');

figure; 
hold on; plot(frcp1{1}.om(1:5:end),frcp1{1}.rho(1:5:end),'kx')
hold on; plot(frcp3{1}.om(1:5:end),frcp3{1}.rho(1:5:end),'rx')

hold on; plot(frcp1{1}.om(1:5:end),frcp1{1}.ampL2(1:5:end),'kx')
hold on; plot(frcp3{1}.om(1:5:end),frcp3{1}.ampL2(1:5:end),'rx')
%
% figure; 
%hold on; plot_stab_lines(pFRCs{2}.om,pFRCs{2}.ampLinf(4,:),pFRCs{2}.st,ST1);
%

%%
% rhosamp = linspace(0,0.5,31); plotdof = [1 2 5];
% zmesh = S.SSM_sensitivity(We,DW,rhosamp,plotdof,pidx,pfracs);