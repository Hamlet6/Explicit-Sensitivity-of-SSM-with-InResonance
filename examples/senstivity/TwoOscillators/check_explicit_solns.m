
demo_explicit;
% also run demo_semi_intrusive to yield result
load('multi_semi.mat') % v2 for first order computation
full(W_0(1).ind)

full(W_0(1).coeffs(:,2))-W10

full(W_0(2).coeffs(:,2))-W11
full(W_0(2).coeffs(:,3))-W20

full(W_0(3).coeffs(:,3))-W21
full(W_0(3).coeffs(:,4))-W30

full(R_0(3).coeffs(1,2))-r21

W_1(2).kappa
W_1(2).W.coeffs-x0
R_1(2).R.coeffs(1)-ftilde
