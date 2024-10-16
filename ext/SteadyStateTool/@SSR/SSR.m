classdef SSR < handle
    %SSR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        M               % Mass matrix
        K               % Stiffness matrix
        C               % Damping matrix (assuming proportional damping)
        T               % Time Period of Oscillation
        Omega           % Frequency basis vector (for quasi-periodic oscillation)
        U               % Matrix of (undamped) eigenvectors
        S               % function handle for the nonlinearity
        DS              % function handle for the Jacobian of the nonlinearity
        f               % function handle for the forcing f(t,T)
        Df              % function handle for derivative of the forcing w.r.t time period T
        order = 1       % order if integration: Newton Cotes is used (existing steps are refined)
        n_steps = 50    % number of intervals in the time domain: 50 by default
        type = 'p'      % 'p': calculations performed for periodic solution (time domain)
                        % 'qp' : calculations performed for periodic solution (time domain)
        %     end
        %
        %     properties (Access = private)
        alpha           % constants defined in (27)
        omega           % constants defined in (27)
        beta            % constants defined in (27)
        gamma           % constants defined in (27)
        omega0          % undamped natural frequencies
        o               % indices of overdamped modes
        c               % indices of critically damped modes
        u               % indices of underdamped modes
        n               % number of DOFs
        t               % time vector 0:dt:T
        dt              % node spacing in the time mesh (T/nt)
        Lvec            % Green's function for displacements (calculated over -T:dt:T)
        DLvec           % Derivative dLdT
        Jvec            % Green's function for velocities (calculated over -T:dt:T)
        ConvMtx         % convolution matrix for displacements
        isupdated       % check if different quantitities are updated according to T
        weights         % Weights vector for higher order integration
        nt              % number of nodes in the time mesh
        Ft              % forcing evaluated over time mesh
        nh              % number of harmonics in of the base frequency Omega
        kappa_set       % the set of frequency-indices used during Fourier transformation
        theta_set       % the set of torus coordinates over which the solution is approximated used during Fourier transformation
        Qmat            % the global Q matrix containing the Green's function (amplification factors) for freq domain analysis
        f_theta         % function handle for forcing in torus coordinates
        f_kappa         % Fourier coefficients for the external forcing
        m               % the discretization along each dimension of the Torus
        E               % The Fourier transform exponential matrix
        Einv            % The inverse Fourier transform exponential matrix
        Evinv           % The inverse Fourier transform exponential matrix
        Ekron           % The Fourier transform exponential matrix kronecker product with Identity
        Einvkron        % The inverseFourier transform exponential matrix kronecker product with Identity
    end
    
    methods
        function O = SSR(M,C,K)
            % Basic system properties
            O.M = M;
            O.C = C;
            O.K = K;
            O.n = length(M);
            
            [V, ~] = eig(full(K),full(M));
            mu = diag(V.' * M * V);
            O.U = V * diag( 1./ sqrt(mu) ); % mass normalized VMs
            O.omega0 = sqrt(diag((O.U.' * K * O.U)));
            zeta = diag((O.U.' * C * O.U))./ (2*O.omega0);
            O.o = find(zeta>1);
            O.c = find(zeta==1);
            O.u = find(zeta<1);
            lambda1 = (-zeta + sqrt(zeta.^2 - 1)).* O.omega0;
            lambda2 = (-zeta - sqrt(zeta.^2 - 1)).* O.omega0;
            O.alpha = real(lambda2);
            O.omega = abs(imag(lambda2));
            O.beta = lambda1;
            O.gamma = lambda2;
        end
        
        function set.T(O,T)
            if ~isequal(O.T,T)
                O.T = T;
                update_t(O);
            end
        end
        
        function set.Omega(O,Omega)
            if ~isrow(Omega)
                Omega = Omega.';
            end
            if ~isequal(O.Omega, Omega)
                O.Omega = Omega;
                not_updated_Q(O);
            end
        end
        
        notupdated(O)
        
        function set.order(O,order)
            if O.order ~= order
                O.order = order;
                update_t(O);
            end
        end
        
        update_t(O)
        
        function set.n_steps(O,n_steps)
            O.n_steps = n_steps;
            update_t(O);
        end        
       
        function set.nh(O,nh)
            if ~isequal(O.nh, nh)
                O.nh = nh;
                not_updated(O);
            end
        end
        
        function set.m(O,m)
            if ~isequal(O.m, m)
                O.m = m;
                not_updated(O);
            end
        end
        
        function theta = get.theta_set(O)
            if ~O.isupdated.theta
                update_theta_set(O);
            end
            theta = O.theta_set;
        end
        
        update_theta_set(O)

        function kappa = get.kappa_set(O)
            if ~O.isupdated.kappa
                update_kappa_set(O);
            end
            kappa = O.kappa_set;
        end
        
        update_kappa_set(O)
        
        function Lvec = get.Lvec(O)
            if ~O.isupdated.L
                update_Lvec(O);
            end
            Lvec = O.Lvec;
        end
        
        function Ft = get.Ft(O)
            if ~O.isupdated.Ft
                update_Ft(O);
            end
            Ft = O.Ft;
        end
        
        update_Ft(O)
        
        function E = get.E(O)
            if ~O.isupdated.E
                update_E(O);
            end
            E = O.E;
        end
        
        update_E(O)
        
        function f_kappa = get.f_kappa(O)
            if ~O.isupdated.f_kappa
                update_f_kappa(O);
            end
            f_kappa = O.f_kappa;
        end
        
        update_f_kappa(O)
        
        update_Lvec(O)
                
        L = L(O,t,T)
        
        function Jvec = get.Jvec(O)
            if O.isupdated.J
                Jvec = O.Jvec;
            else
                update_Jvec(O);
                Jvec = O.Jvec;
            end
        end
        
        update_Jvec(O)
        
        J = J(O,t,T)
        
        CM = get_ConvMtx(O)
        
        update_ConvMtx(O)
        
        update_dLvec(O)
        
        L = dLdT(O,t,T)
        
        function Q = get.Qmat(O)
            if ~O.isupdated.Q
                update_Qmat(O);
            end
            Q = O.Qmat;
        end
        
        update_Qmat(O)
        
        eta = convolution_x(O,phi)
        
        etad = convolution_xd(O,phi)
        
        eta = convolution_dT(O,phi)
        
        F_P = F_P(O,x)
        
        F_Q = F_Q(O,x)
        
        DF_P = DF_P(O,x)
        
        DF_Q = DF_Q(O,x)
        
        [x,xd] = LinearResponse(O)
        
        [x, xd] = Picard(O,varargin)
        
        [x, xd] = NewtonRaphson(O,varargin)
        
        [x0,tol,maxiter] = parse_iteration_inputs(O,inputs)
        
        [OMEGA, SOL, PICARD] = sequential_continuation(O,range,varargin)
    end
    
    methods(Static)
        w = NewtonCotes(order)
        
        h = h(t)
        
        F = evaluate_fun_over_array(f,x,vectorized)
    end
end

