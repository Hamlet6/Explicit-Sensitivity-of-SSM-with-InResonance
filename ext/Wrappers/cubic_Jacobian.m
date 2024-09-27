function f3 = cubic_Jacobian(Assembly,x,y,z)
% this function computes the Jacobian of the cubic component of the
% nonlinear internal force in global coordinates at the element
% level. The Jacobian is evaluated along the direction (x,y)
% and acted on the vector z. Hence, the output is a vector.
X = Assembly.unconstrain_vector([x y z]);
f3 = Assembly.constrain_vector(Assembly.vector('DF3',X(:,1),X(:,2),X(:,3)));
end