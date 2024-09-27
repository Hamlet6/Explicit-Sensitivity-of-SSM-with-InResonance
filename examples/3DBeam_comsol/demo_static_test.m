
clear all
% close all
clc
%% generate model

model = mphload('model_3d.mph');
out = mphgetu(model,'soltag','sol1');
%%

ML = mphmatrix(model ,'sol1','Out', {'L','K','ud','Nullf'});
L = ML.Nullf'*(ML.L-ML.K*ML.ud);
outdof = find(L);

model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);

MA = mphmatrix(model ,'sol1', ...
 'Out', {'K','L','D','E','Null','Nullf','ud','uscale'});
M_0 = MA.E; C_0 = MA.D; K_0 = MA.K; L_0 = MA.L; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud; uscale = MA.uscale;

M = Nullf'*M_0*Null;
C = Nullf'*C_0*Null;
K = Nullf'*K_0*Null;

n = length(M);
disp(['Number of degrees of freedom = ' num2str(n)])
disp(['Phase space dimensionality = ' num2str(2*n)])
[fint,dfint] = get_fint(K,model,Null, Nullf, ud);
%%
x0 = zeros(n,1);

load = L;
func = @(x) K*x + fint(x) - load;
func_linear = @(x) K*x - load;
% u + ud;
u = Null*newton_raph(func, x0, dfint, K) + ud;
u_linear = Null*newton_raph(func_linear, x0, dfint, K) + ud;

diff1_norm = norm( out - u )/norm(out)
diff2_norm = norm( u - u_linear )/norm(out)

%%
% yn = newton_raph(func, x0, dfint, K);
% yl = newton_raph(func_linear, x0, dfint, K);

%% plot deflection

u1 = out(3:60:end);
u2 = u(3:60:end);
u3 = u_linear(3:60:end);
x = linspace(0,1,numel(u1));
figure
plot(x,u1,'DisplayName','COMSOL solution','linewidth',2)
hold on
plot(x,u2,'--','DisplayName','Nonlinear','linewidth',2)
hold on
plot(x,u3,'-o','DisplayName','Linear','linewidth',2)
legend
xlabel('scaled length x/L')
ylabel('Displacement D')
%% linear analysis results
Len = 2; h=0.1; E = 71e9; nv = 0.33; b = 0.16;
I = b*h^3/12;
d = norm(load,'inf')*Len^3/(3*E*I)
norm(u,'inf')

