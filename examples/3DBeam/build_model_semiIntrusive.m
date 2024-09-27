function [M,C,K,fnl, fext, outdof,MyAssembly] = build_model_semiIntrusive(elementType,nx,ny,nz)
run ../../install.m
%% PREPARE MODEL
%% Read Mesh
startLIN = tic;

%% PREPARE MODEL

% DATA ____________________________________________________________________
E       = 70e9;     % Young's modulus [Pa]
rho     = 2700;     % density [kg/m^3]
nu      = 0.33;     % Poisson's ratio

% MESH_____________________________________________________________________
l = 3;
w = .3;
t = .1;

[nodes, elements, nset]=mesh_3Dparallelepiped(elementType,l,w,t,nx,ny,nz);

% % Alternatively, one can write an input file in ABAQUS and read it as:
% filename = 'Job-BeamHex';
% [nodes, elements, nset, elset] = mesh_ABAQUSread(filename); % HEX20 mesh
PlotMesh(nodes,elements)
%% FE model
disp('Building FE model')
% Material
myMaterial = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
% Element
switch elementType
    case 'HEX20'
        myElementConstructor = @()Hex20Element(myMaterial);
    case 'TET10'
        myElementConstructor = @()Tet10Element(myMaterial);
end

MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

%% Assemble linear stiffness, mass and damping
disp('Assembling M,C,K matrices')
% % parallelized assembly

MyAssembly = Assembly(MyMesh);

u0 = zeros(MyMesh.nDOFs, 1);

filename = ['matrices_' num2str(MyMesh.nElements) '.mat'];


[K,~] = MyAssembly.tangent_stiffness_and_force(u0);
M = MyAssembly.mass_matrix();
% C = MyAssembly.damping_matrix();
MyAssembly.DATA.K = K;
MyAssembly.DATA.M = M;
%% apply boundary conditions
disp('Applying boundary conditions')

% MESH > BOUNDARY CONDITONS
MyMesh.set_essential_boundary_condition([nset{1} nset{4}],1:3,0)
% MyMesh.BC.set_dirichlet_dofs([nset{2} nset{3}],1:3,0) % abaqus

M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);
% C = MyAssembly.constrain_matrix(C);


%% Eigenvalue problem
disp('Solving undamped eigenvalue problem')
n_VMs = 3; % first n_VMs modes with lowest frequency calculated
[V0,omega2] = eigs(K,M,n_VMs,'SM');
omega = sqrt(diag(omega2));


% [f0,ind] = sort(sqrt(diag(omega2))/2/pi);
% V0 = V0(:,ind);
% for ii = 1:n_VMs
%     V0(:,ii) = V0(:,ii)/max(sqrt(sum(V0(:,ii).^2,2)));
% end

V = MyAssembly.unconstrain_vector(V0);
mod = 1;
v1 = reshape(V(:,mod),3,[]).';

%{
figure('units','normalized','position',[.2 .1 .6 .8])
h = PlotFieldonDeformedMesh(nodes,elements,v1,'factor',10);
set(h,'Linestyle','none')
colormap jet
title(['Mode ' num2str(mod) ', Frequency = ' num2str(omega(mod)/(2*pi),3) ' Hz'] )
drawnow
light
%}

%% Damping matrix
disp('Using Rayleigh damping')
W =   omega(1:2);
a = [W(1) 1/W(1);W(2) 1/W(2)]\[0.004;0.004];
C = a(2) * M + a(1) * K;

MyAssembly.DATA.C = a(2) * MyAssembly.DATA.M + a(1) * MyAssembly.DATA.K;    fnl = cell(1,2);
computationTimeLIN = toc(startLIN);

save(filename,'M','C','K','computationTimeLIN','-v7.3')

%% external force assembly
disp('Assembling external force vector')

% nForceElements = numel(nbcElements);
% weights = sparse(nbcElements,ones(nForceElements,1),true,MyMesh.nElements,1);
% fext = MyAssembly.uniform_body_force(weights);
% fext = MyAssembly.constrain_vector(fext);


outcoord = [l/2,w/2,t/2];
dist = zeros(MyMesh.nNodes,1);

for j = 1:MyMesh.nNodes
    dist(j) = norm(nodes(j,1:3) - outcoord);
end

[~,outnode] = min(dist);
outdof = outnode*3; % x direction

% forcing tip nodes
fext = 10^5*sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
fext = MyAssembly.constrain_vector(fext);
outdof = find(fext);

%% Tensor Assembly
%% nonlinearity
disp('Getting nonlinearity coefficients')

fnl{1} = @(input)quadratic_nonlinearity(MyAssembly,input{:});
fnl{2} = @(input)cubic_nonlinearity(MyAssembly,input{:});
%
%dfnl{1} = @(input)quadratic_Jacobian(MyAssembly,input{:});
%dfnl{2} = @(input)cubic_Jacobian(MyAssembly,input{:});


end

