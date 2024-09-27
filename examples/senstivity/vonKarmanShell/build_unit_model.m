function [M,K,fnl,fext,outdof] = build_unit_model(nDiscretization)


%% Finite Element Setup
% Geometry
l = 1; % length of domain [m]
b = 2; % breadth of domain [m]
t = 1e-2; % thickness of plate [m]
w = 1e-1; % curvature parameter (height of the midpoint relative to ends) [m]
% material properties
E       = 1;  % 70e9 % 200e9 % Young's modulus [Pa]
rho     = 1;  % 2700 % 7850 % density [kg/m^3]
nu      = 0.33;    % Poisson's ratio 
kappa   = 1e5; % material damping modulus 1e8

%% FE model
disp('Building FE model')
% Material
myMaterial  = KirchoffMaterial();
set(myMaterial,'YOUNGS_MODULUS',E,'DENSITY',rho,'POISSONS_RATIO',nu,'DAMPING_MODULUS',kappa);
% Element
myElementConstructor =  @()TriShellElement(t, myMaterial); % same element all across the domain

% Meshing the geometry
nl = nDiscretization;
nb = 2*nDiscretization; 
[nodes,elements,bnodes] = RectangularMesh(l,b,nl,nb,w);     % Rectangular Mesh definition

% creating Mesh object
MyMesh = Mesh(nodes);
MyMesh.create_elements_table(elements,myElementConstructor);

% Plot mesh
figure('Name','Mesh'); PlotMesh(nodes,elements,0);

disp('Getting nonlinearity coefficients')
filename = ['tensors_' num2str(MyMesh.nElements) '.mat'];
try     
    load(filename,'M','K','fnl','outdof','fext');
    disp('Loaded matrices and tensors from storage');
catch    
    %% Assemble linear stiffness, mass and damping
    
    MyAssembly = Assembly(MyMesh);
    K = MyAssembly.stiffness_matrix();
    M = MyAssembly.mass_matrix();
    
    %% apply boundary conditions
    disp('Applying boundary conditions')
    MyMesh.set_essential_boundary_condition([bnodes{3}, bnodes{4}],1:3,0) % simply supported on opposite ends
    M = MyAssembly.constrain_matrix(M);
    K = MyAssembly.constrain_matrix(K);    
    
    %% Eigenvalue problem
    disp('Solving undamped eigenvalue problem')
    n_VMs = 5; % first n_VMs modes with lowest frequency calculated 
    [V0,omega2] = eigs(K,M,n_VMs,'SM');
    omega = sqrt(diag(omega2));
    
    V = MyAssembly.unconstrain_vector(V0);
    mod = 1;
    v1 = reshape(V(:,mod),6,[]);
    figure;
    PlotFieldonDeformedMesh(nodes,elements,v1(1:3,:).','factor',5);
    title(['Frequency = ' num2str(omega(mod)/(2*pi)) ' Hz'] )
    set(colorbar,'visible','off')
    
    %% external force assembly
    disp('Assembling external force vector')
    outcoord = [l/2,b/4]; % output coordinate
    outdir = 3; % transverse displacement
    dist = vecnorm(MyMesh.nodes(:,1:2) - repmat(outcoord,[MyMesh.nNodes,1]),2,2);
    [~,outnode] = min(dist);
    outdof = (outnode-1)*MyAssembly.Mesh.nDOFPerNode+outdir;
    
    outdofvec = sparse(outdof,ones(size(outdof)),1,MyMesh.nDOFs,1);
    outdofvec = MyAssembly.constrain_vector(outdofvec);
    outdof = find(outdofvec);
    
    fext = 100*outdofvec;    
    %% Tensor Assembly
    fnl = cell(1,2);
    disp('Assembling Tensors')
    fnl{1} = MyAssembly.tensor('T2',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3]);
    fnl{2} = MyAssembly.tensor('T3',[MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs, MyMesh.nDOFs], [2,3,4]);           
    % apply boundary conditions
    for j = 1:length(fnl)
        fnl{j} = MyAssembly.constrain_tensor(fnl{j});
    end  
    save(filename,'M','K','fnl','fext','outdof');
end
end
