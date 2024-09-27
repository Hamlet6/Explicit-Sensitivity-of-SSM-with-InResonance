function fint = internalForce(Assembly,x)
% this function computes the Jacobian of the cubic component of the
% nonlinear internal force in global coordinates at the element
% level. 

x = Assembly.unconstrain_vector(x);
fint = Assembly.constrain_vector(Assembly.vector('internal_force',x));
end