function f2 = quadratic_Jacobian(Assembly,x,y)
% this function computes the Jacobian of the quadratic component
% of the nonlinear internal force in global coordinates
% at the element level. The Jacobian is evaluated along the
% direction (x) and acted on the vector y. Hence, the output is a vector.

x = Assembly.unconstrain_vector(x);
y = Assembly.unconstrain_vector(y);
f2 = Assembly.constrain_vector(Assembly.vector('DF2',x,y));
end