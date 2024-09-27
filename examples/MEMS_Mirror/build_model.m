function [M,C,K,fnl,fext,outdof] = build_model(forcing)


%% Read Mesh
startLIN = tic;
mshfile = 'Mirror.msh2';
disp(['Reading mesh from ' mshfile]) 
[nodes, elements, constrainedNodes, nbcElements] = read_msh(mshfile);

fprintf('Object dimension is %d x %d x %d um \n' , max(2*abs(nodes(:,1))),max(2*abs(nodes(:,2))),max(abs(nodes(:,3))))

%% Rescaling nodes to micrometers

% centering z axis - so center of structure is at (0,0,0)

thickness = max(nodes(:,3));
nodes(:,3) = nodes(:,3) - thickness/2;

scale = 1e-6; % units of model given in micrometers
nodes = nodes*scale;


%% parameters

% MEMS Gyro - Silicon
E       = 148e9;   	% Young's modulus [GPa] (in m)
rho     = 2.3e3; 	% density [kg/m^3]
nu      = 0.23;     % Poisson's ratio

%% FE model
disp('Building FE model')

% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu);
% Element
myElementConstructor = @()Hex8Element(myMaterial);

% Mesh object
MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

% Plot mesh
figure('Name','Mesh'); PlotMesh(nodes,elements,0);

%% Assemble linear stiffness, mass and damping
disp('Assembling M,C,K matrices')
% % parallelized assembly
% cluster = parcluster('local');
% cluster.NumWorkers = 4;
% parpool(cluster, 4)
% MyAssembly = Assembly(myMesh,true); 
MyAssembly = Assembly(MyMesh);



u0    = zeros(MyMesh.nDOFs, 1);
[K,~] = MyAssembly.tangent_stiffness_and_force(u0);
M     = MyAssembly.mass_matrix();

%% apply boundary conditions
disp('Applying boundary conditions')
MyMesh.set_essential_boundary_condition(constrainedNodes,1:3,0) % fixed on the ends of torsionbars
%}
M = MyAssembly.constrain_matrix(M);
K = MyAssembly.constrain_matrix(K);

%% Eigenvalue problem

disp('Solving undamped eigenvalue problem')
n_VMs = 10; % first n_VMs modes with lowest frequency calculated 
[V0,omega2] = eigs(K,M,n_VMs,'sm');
omega = sqrt(diag(omega2));

% Plot modes
mod = 1;

V = MyAssembly.unconstrain_vector(V0);
v1 = reshape(V(:,mod),3,[]).';

figure('units','normalized','position',[.2 .1 .6 .8])
h = PlotFieldonDeformedMesh(nodes,elements,v1,'factor',1e-2);
set(h,'Linestyle','none')
colormap jet
title(['Mode ' num2str(mod) ', Frequency in rad = ' num2str(omega(mod))] )
drawnow
light


%% Damping matrix
disp('Using Rayleigh damping')
W =   omega(1:2);

a =  [W(1) 1/W(1);W(2) 1/W(2)]\[1;1]*1e-5

C = a(2) * M + a(1) * K;
%% external force assembly

disp('Assembling external force vector')

% Torsion bar along y axis
% Mirror along y axis. Outdof node at x ~ 0 and outmost in y direction
outnode = 316; 
forcenode = [316,317]; % adjacent nodes, force symmetric with respect to torsion bar axis  

fprintf('Forcing the nodes at')
nodes(forcenode,:)

fprintf('################################################### \n')

switch forcing
    case 'x'
        fprintf('Forcing the x coodinate!')
        forcingdofs = (forcenode-1)*MyAssembly.Mesh.nDOFPerNode + 1;
        outdof = (outnode-1)*MyAssembly.Mesh.nDOFPerNode+1;

    case 'y'
        fprintf('Forcing the y coodinate!')
        forcingdofs = (forcenode-1)*MyAssembly.Mesh.nDOFPerNode + 2;
        outdof = (outnode-1)*MyAssembly.Mesh.nDOFPerNode+2;

    case 'z'
        fprintf('Forcing the z coodinate!')
        forcingdofs = (forcenode-1)*MyAssembly.Mesh.nDOFPerNode + 3;
        outdof = (outnode-1)*MyAssembly.Mesh.nDOFPerNode+3;

end

% get outdof
% forcing tip nodes
outdof = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
outdof = MyAssembly.constrain_vector(outdof);
outdof = find(outdof);


% forcing tip nodes
fext = sparse(forcingdofs,ones(size(forcingdofs)),1,MyMesh.nDOFs,1);
fext = MyAssembly.constrain_vector(fext);

computationTimeLIN = toc(startLIN)

%% Tensor Assembly


disp('Getting nonlinearity coefficients')
filename = ['tensors_' num2str(MyMesh.nElements) '.mat']
try     
    load(filename,'fnl')
    disp('Loaded tensors from storage')
    load(filename, 'computationTimeTensors')
    disp(['Total time spent on model assembly = ' datestr(datenum(0,0,0,0,0,computationTimeTensors + computationTimeLIN),'HH:MM:SS')])
catch
    fnl = cell(1,2);
    disp('Assembling Tensors')
    startTensors = tic;
    fnl{1} = MyAssembly.tensor('T2',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3]);
    fnl{2} = MyAssembly.tensor('T3',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3,4]);       
    computationTimeTensors = toc(startTensors);

    fprintf('Tensors are constrained before saving')
    for j = 1:length(fnl)
        fnl{j} = MyAssembly.constrain_tensor(fnl{j});
    end
    disp('Saving Tensors')
    save(filename,'fnl','computationTimeTensors','-v7.3')
    disp(['Total time spent on model assembly = ' datestr(datenum(0,0,0,0,0,computationTimeTensors + computationTimeLIN),'HH:MM:SS')])
end
