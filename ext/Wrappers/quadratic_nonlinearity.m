function f2 = quadratic_nonlinearity(Assembly,x,y)
x = Assembly.unconstrain_vector(x);
y = Assembly.unconstrain_vector(y);
f2 = Assembly.constrain_vector(Assembly.vector('F2',x,y));
end