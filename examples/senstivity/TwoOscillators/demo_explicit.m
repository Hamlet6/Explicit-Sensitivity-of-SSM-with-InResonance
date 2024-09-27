%% 
% Consider the following system with 2 dofs
% 
% $$\ddot{x}_1+c_1\dot{x}_1+x_1+b_1x_1x_2+x_1^3+x_1x_2^2=\epsilon f_1\cos\Omega 
% t\\\ddot{x}_2+c_2\dot{x}_2+2x_2+b_2x_1^2+0.3x_1x_2^2+1.6x_2^3=\epsilon f_2\cos\Omega 
% t$$

clear all
%% Setup model

m = 1;
c1 = 5e-3;
c2 = 1e-2;
% c1 = 0.02; c2 = 0.02;
b1 = 1;
b2 = 1;
f1 = 1;
f2 = 0;
[mass,damp,stiff,fnl,fext] = build_model_semi(c1,c2,b1,b2,f1,f2);

DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl_semi',fnl);
set(DS.Options,'Emax',1,'Nmax',2)


[V,D,W] = DS.linear_spectral_analysis();
%% SSM computation

S = SSM(DS);
resonant_modes = [1 2]; % choose master spectral subspace
S.choose_E(resonant_modes);
order = 3;                  % SSM expansion order

set(DS,'fext',fext,'Omega',0.98*imag(D(1)));
%% explicit computation
tic
[W10,W20,W11,W30,W21,r21,x0,ftilde] = S.explicit_whisker();
toc


%% FIRST ORDER COMPUTATION
DS = DynamicalSystem();
% set(DS,'B',[damp mass;mass zeros(2)],'A',[-stiff zeros(2);zeros(2) mass],'fnl_semi',fnl);
set(DS,'M',mass,'C',damp,'K',stiff,'fnl_semi',fnl, 'order', 1);
set(DS.Options,'Emax',5, 'Nmax',10)

[V,D,W] = DS.linear_spectral_analysis();
%% SSM computation

S = SSM(DS);
resonant_modes = [1 2]; % choose master spectral subspace
S.choose_E(resonant_modes);
order = 3;                  % SSM expansion order

%% explicit computation
tic
[W10,W20,W11,W30,W21,r21] = S.explicit_whisker();
toc

