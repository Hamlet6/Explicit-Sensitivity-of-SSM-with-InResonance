function dfint = internalForceJacobian(Assembly,x)
% this function computes the Jacobian of the cubic component of the
% nonlinear internal force in global coordinates at the element
% level. 

x = Assembly.unconstrain_vector(x);
% only takes first output argument which is the tangent stiffness
% this also includes the linear system matrix K
[dfint] = Assembly.constrain_matrix(Assembly.matrix('tangent_stiffness_and_force',x));
end