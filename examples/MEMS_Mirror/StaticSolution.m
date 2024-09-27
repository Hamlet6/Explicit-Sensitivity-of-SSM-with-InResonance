%% MYMEMS
clearvars
clc
close all

%% Read Mesh
mshfile = 'Mirror.msh2';
disp(['Reading mesh from ' mshfile]) 
[nodes, elements, constrainedNodes, nbcElements] = read_msh(mshfile);

% due to torsion bar, the mirror is slightly wider than 600 um.
% due to the volumes which are used for boundary condition at end of
%         torsionbar the x-dimension is 2x 0.1um longer than 1400 um.

fprintf('Object dimension is %d x %d x %d um \n' , max(2*abs(nodes(:,1))),max(2*abs(nodes(:,2))),max(abs(nodes(:,3))))

% Rescale if calculation is carried out in meters
scale = 1e-6; % units of model given in micrometers
nodes = nodes*scale;

% set center of structure to (0,0,0)

nodes(:,3) = nodes(:,3) - max( nodes(:,3))/2;
%% parameters
% Silicon
E       = 148e9;   	% Young's modulus [GPa] (in m)
rho     = 2.3e3; 	% density [kg/m^3]
nu      = 0.23;     % Poisson's ratio

%% FE model
disp('Building FE model')
% Material
myMaterial           = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
myElementConstructor = @()Hex8Element(myMaterial);

% Mesh object
MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

%% Assemble linear stiffness, mass and damping
disp('Assembling K matrix') 
MyAssembly  = Assembly(MyMesh);
u0          = zeros(MyMesh.nDOFs, 1);
[K,~]       = MyAssembly.tangent_stiffness_and_force(u0);

% Boundary conditions
MyMesh.set_essential_boundary_condition(constrainedNodes,1:3,0) % cantilever on one edge
K           = MyAssembly.constrain_matrix(K);

% Plot mesh
figure('Name','Mesh'); PlotMesh(nodes,elements,0);

%% Get Forcing Nodes

% Torsion bar along y axis
% Mirror along y axis.

% All nodes on top of positive y side of mirror
pos_y_idx = find(max(nodes(:,2)) == nodes(:,2));
pos_idx   =pos_y_idx( find(max(nodes(pos_y_idx,3)) == nodes(pos_y_idx,3)));

% All nodes on bottom of negative y side of mirror
neg_y_idx = find(min(nodes(:,2)) == nodes(:,2));
neg_idx   = neg_y_idx( find(min(nodes(neg_y_idx,3)) == nodes(neg_y_idx,3)) );

fprintf('Nodes with positive forcing')
pos_nodes = nodes(pos_idx,:);
fprintf('Nodes with negative forcing')
neg_nodes = nodes(neg_idx,:);


% Convert nodes to outdof
fdof  =  3* ([pos_idx;neg_idx] -1)  + 3; % z direction
fdofpos =  3* ([pos_idx] -1)  + 3; % z direction
fdofnneg =  3* ([neg_idx] -1)  + 3; % z direction

%% Vector of meshpoints

% get forcing 
fext = 1e-6 * sparse(fdof,ones(size(fdof)),[ones(size(fdof,1)/2,1),-ones(size(fdof,1)/2,1)],MyMesh.nDOFs,1);
fext = MyAssembly.constrain_vector(fext);

% get constrained coordinates of outdof
fdofpos_restr = sparse(fdofpos,ones(size(fdofpos)),1,MyMesh.nDOFs,1);
fdofpos_restr = MyAssembly.constrain_vector(fdofpos_restr);
fdofpos_restr = find(fdofpos_restr);

fdofneg_restr =sparse(fdofnneg,ones(size(fdofnneg)),1,MyMesh.nDOFs,1);
fdofneg_restr = MyAssembly.constrain_vector(fdofneg_restr);
fdofneg_restr = find(fdofneg_restr);

%% Solve static solution

% Static solution to antisymmetric forcing
static = K \ fext ;

fprintf('Static solution to antisymmetric forcing \n')
fprintf('Static displacement at negative y end of mirror')
static(fdofneg_restr)
fprintf('Static displacement at positive y end of mirror')
static(fdofpos_restr)
fprintf('Difference of the absolute value of opposite displacements')
full([static(fdofpos_restr) + static(fdofneg_restr)])
fprintf('relative difference')
full([(static(fdofpos_restr) + static(fdofneg_restr))./static(fdofneg_restr)])


V = MyAssembly.unconstrain_vector(static);
v1 = reshape(V,3,[]).';

figure()
h = PlotFieldonDeformedMesh(nodes,elements,v1,'factor',1);
colormap jet
title(['Static solution '] )
drawnow