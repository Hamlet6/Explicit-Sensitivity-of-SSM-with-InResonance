function [M,C,K,fnl,dfnl2,dfnl3,outdof] = build_model_semi(nElements,kLin,kNon)
%% Finite Element Setup
assert(mod(nElements,2)==0, 'The number of elements is not an even number.');
% Geometry

h=10;                                   %Height of beam [mm]
b=10;                                   %Width of beam [mm]
l=2700;                                 %Length of beam [mm]
% Mesh parameters

% Material properties
E=45000000;                             %Young's modulus [kPa]   
rho=1780*10^(-9);                       %Density [kg/mm^3]    
kappa = 1e3;
nu    = 0.3;    % nu

%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
% Element
myElementConstructor = @()BeamElement(b, h, myMaterial); % same element all across the domain

% Meshing the geometry
dx = l/nElements;
x = (0:dx:l).';
nNodes = size(x,1);
nodes = [x, zeros(nNodes,1)];
elements = [1:nNodes-1;2:nNodes].';

% creating Mesh object
MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

%% Assemble linear stiffness, mass and damping
disp('Assembling M,C,K matrices')
MyAssembly = Assembly(MyMesh);
K = MyAssembly.stiffness_matrix();
M = MyAssembly.mass_matrix();
C = MyAssembly.damping_matrix();

%% apply boundary conditions
disp('Applying boundary conditions')
MyMesh.set_essential_boundary_condition(1,[1 2 3],0);     % clamped at left end
MyMesh.set_essential_boundary_condition(nNodes, [1 2],0); % hinged at right end
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
C = MyAssembly.constrain_matrix(C);

%% external force assembly
disp('Assembling external force vector')

outnode = [1 2]/4*nElements+1;
outdof = outnode(2)*3-1; % transverse direction

outdofvec = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
outdofvec = MyAssembly.constrain_vector(outdofvec);
outdof = find(outdofvec);

K(outdof,outdof) = K(outdof,outdof)+kLin;

%% nonlinearity
disp('Getting nonlinearity coefficients')

fnl{1} = @(input) quadratic_nonlinearity(MyAssembly,input{:});
fnl{2} = @(input) cubic_nonlinearity(MyAssembly,input{:})+cubic_support(input,kNon,outdof);
n     = size(M,1);
dfnl2 = @(input) zeros(n,1);
dfnl3 = @(input) dcubic_support(input,outdof); 
dfnl2 = {dfnl2, dfnl2}; % w.r.t kLin and kNon
dfnl3 = {@(input) zeros(n,1), dfnl3};

end


function y = cubic_support(x,kNon,outdof)
v1 = x{1};
v2 = x{2};
v3 = x{3};
y  = zeros(size(v1));
y(outdof) = kNon*v1(outdof)*v2(outdof)*v3(outdof);
end

function y = dcubic_support(x,outdof)
v1 = x{1};
v2 = x{2};
v3 = x{3};
y  = zeros(size(v1));
y(outdof) = v1(outdof)*v2(outdof)*v3(outdof);
end
