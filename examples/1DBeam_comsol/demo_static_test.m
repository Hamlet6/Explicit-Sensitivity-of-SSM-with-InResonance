
clear all
% close all
clc
%% generate model

model = mphload('model_1d.mph');
out = mphgetu(model,'soltag','sol1');

%%

ML = mphmatrix(model ,'sol1','Out', {'L','K','ud','Nullf'});
L = ML.Nullf'*(ML.L-ML.K*ML.ud);
outdof = find(L); 
% we deduce from outdof that the ordering is flipped, namely the dofs at
% the free end is put at the first place. In addition, we know that the
% dofs at each node is in the ordering (thx,thy,thz,ux,uy,uz).

model.component('comp1').physics('beam').feature('pl1').set('Fp', [0 0 0]);

MA = mphmatrix(model ,'sol1', ...
 'Out', {'K','Kc','L','D','E','Null','Nullf','ud','uscale'});
M_0 = MA.E; C_0 = MA.D; K_0 = MA.K; L_0 = MA.L; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud; uscale = MA.uscale;

% Null  = diag(uscale)*Null;
% Nullf = diag(uscale)*Nullf;
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
u1 = out(5:6:end);
u2 = u(5:6:end);
u3 = u_linear(5:6:end);
x = linspace(1,0,numel(u1));
figure
plot(x,u1,'DisplayName','COMSOL solution','linewidth',2)
hold on
plot(x,u2,'--','DisplayName','Nonlinear','linewidth',2)
hold on
plot(x,u3,'-o','DisplayName','Linear','linewidth',2)
legend
xlabel('scaled length x/L')
ylabel('Displacement D')
%%
figure;
plot(x,sqrt(out(4:6:end).^2+out(5:6:end).^2));
%% linear analysis results
Len = 1; h=0.1; E = 71e9; nv = 0.33; b = 0.16;
I = b*h^3/12;
d = norm(load,'inf')*Len^3/(3*E*I)
th = norm(load,'inf')*Len^2/(2*E*I)
