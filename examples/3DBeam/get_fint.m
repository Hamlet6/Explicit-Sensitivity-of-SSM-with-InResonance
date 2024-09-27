function [fnl] = get_fint(K,elementType,nx,ny,nz)

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio

% MESH_____________________________________________________________________
l = 3;
w = .3;
t = .1;

[nodes, elements, nset]=mesh_3Dparallelepiped(elementType,l,w,t,nx,ny,nz);

PlotMesh(nodes,elements)
%% FE model
disp('Building FE model')
% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
% Element

myElementConstructor = @()Tet10Element(myMaterial);


MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);


MyAssembly = Assembly(MyMesh);
MyMesh.set_essential_boundary_condition([nset{1} nset{4}],1:3,0)


fnl = @(input) internalForce(MyAssembly,input) - K*input; % 


end

