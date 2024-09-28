%% In this demo, we check the equivalence of function handle evaluation and 
% tensor evaluation

clear all; 
% system parameters
m = 1;
c1 = 5e-3;
c2 = 1e-2;
b1 = 1;
b2 = 1;
f1 = 1;
f2 = 0;

% generate model
[mass,damp,stiff,fnl,fext] = build_model_semi(c1,c2,b1,b2,f1,f2);
[mass,damp,stiff,fnl2,fext] = build_model(c1,c2,b1,b2,f1,f2);

% quadratic
v1 = rand(numel(fext),1);
v2 = rand(numel(fext),1);

quad_hand = fnl{1}({v1;v2});  % by function handle
quad_tens = double(ttv(fnl2{1},{v1,v2},[2,3])); % by tensor operation
quad_err  = norm(quad_tens-quad_hand)/norm(quad_tens); % residual

fprintf('relative error at quadratic order is %d\n',quad_err);

% cubic
v3 = rand(numel(fext),1);

cubic_hand = fnl{2}({v1;v2;v3});  % by function handle
cubic_tens = double(ttv(fnl2{2},{v1,v2,v3},[2,3,4])); % by tensor operation
cubic_err  = norm(cubic_tens-cubic_hand)/norm(cubic_tens); % residual

fprintf('relative error at cubic order is %d\n',cubic_err);
