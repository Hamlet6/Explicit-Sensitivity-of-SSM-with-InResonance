L = 1; h = 0.1; E = 71e9; nv = 0.33; b = 0.16;
I = b*h^3/12;
F = 20000;
d = F*L^3/(3*E*I)

% cantilever beam frequency
rho = 2.71e3; 
omega = [1.875; 4.694; 7.855].^2*sqrt(E*I/(rho*b*h*L^4))

Is = h*b^3/12; sqrt(Is/I)
omega = [1.875; 4.694; 7.855].^2*sqrt(E*Is/(rho*b*h*L^4))


E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio

% MESH_____________________________________________________________________
l = 3;
w = .3;
t = .1;

F = epsilon*norm(f_0,'inf');
I = w*t^3/12;
delta = F*l^3/(192*E*I)

