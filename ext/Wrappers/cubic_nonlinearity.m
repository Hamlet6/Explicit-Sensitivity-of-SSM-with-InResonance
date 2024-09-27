function f3 = cubic_nonlinearity(Assembly,x,y,z)
X = Assembly.unconstrain_vector([x y z]);
f3 = Assembly.constrain_vector(Assembly.vector('F3',X(:,1),X(:,2),X(:,3)));
end