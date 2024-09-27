
nDiscretization = 10;
E = 1; rho = 1; % unit value
[M1,C1,K1,fnl1,fext1,outdof1] = build_model(nDiscretization,E,rho,'type1');

E = 123; rho = 157;
[M2,C2,K2,fnl2,fext2,outdof2] = build_model(nDiscretization,E,rho,'type1');

% check proportionality relation
res1 = M2-M1*rho; norm(res1(:),'inf')
res2 = K2-K1*E;   norm(res2(:),'inf')

norm(fnl2{1}-fnl1{1}*E)/norm(fnl2{1})
norm(fnl2{2}-fnl1{2}*E)/norm(fnl2{2})