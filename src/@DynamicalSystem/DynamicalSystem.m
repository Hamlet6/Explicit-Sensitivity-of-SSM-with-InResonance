classdef DynamicalSystem < matlab.mixin.SetGetExactNames
    % DynamicalSystem construct a dynamical system object in first order or
    % second order form
    
    % M\ddot{x} + C\dot{x} + K x + fnl(x,xd) = fext(t) - Second order
    % B \dot{z} = F(z) + Fext(t)                    - First order
    % Here fnl(x) is a polynomial function of degree two or higher, which
    % is stored as a cell array such that fnl{k} corresponds to polynomials
    % of degree k+1. fnl{k} is given by a tensor of order k+2, where the
    % first mode corresponds to indices for the force vector.
    % Likewise, F(z) is a polynomial function of degree one or higher,
    % i.e., F(z) = Az + Higher order terms. F is stored as a cell array,
    % where the i-th entry gives the tensor/multiindex representation of
    % polynomials of degree i.
    
    % The second order form is converted into the first order form with z =
    % [x;\dot{x}], B = [C M;M 0], A = [-K 0;0 M], F(z)=Az+[-fnl(x,xd);0], and
    % Fext(t) = [fext(t); 0]
    
    properties
        M = []
        C = []
        K = []
        A = []
        B = []
        BinvA
        % senstivity info
        DM = [] % derivative w.r.t parameters
        DK = [] % derivative w.r.t parameters
        DC = [] % derivative w.r.t parameters
        % Rayleigh damping info C = alpha(om)*M+beta(om)*K
        al  = [] 
        be  = []
        dal = [] % derivative w.r.t omega
        dbe = [] % derivative w.r.t omega
        
        %% Internal Forces
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Intrusive Nonlinearity
        fnl = []            % second order nonlinear internal forces
        F   = []            % first order nonlinear internal forces
                
        %   Semi-intrusive Nonlinearity        
        fnl_semi    = []    % function handles and parameters to compute nonlinearity, has to take N dim vectors as input
        dfnl_semi   = []    % Handle for nonlinearity derivative
        Dfnl2_semi  = []    % derivative of quadratic function w.r.t parameters
        Dfnl3_semi  = []    % derivative of cubic function w.r.t parameters
                
        F_semi  = []        % function handles for first order dynamical system
        dF_semi = []        % function handle for nonlinearity derivative
        F_semi_sym  = true  % by default a symmetric function handle is assumed        
        
        %   Non-intrusive Nonlinearity
        fnl_non = [];       % second order internal nonlinear force function                     
        dfnl_non  = [];     % second order internal nonlinear force jacobian function
        
        F_non   = [];       % first order internal nonlinear force function
        dF_non  = [];       % first order internal nonlinear force jacobian function
        
        nl_damp = false     % indicates whether 2nd order fct handles for nonlinearities need
                            % only displacement or full phase space vectors as input
        %% Other parameters
        fext = []
        Fext = []
        Omega = []
        dfext = []          % derivative of fext w.r.t parameters
        
        n                   % dimension for x
        N                   % dimension for z
        order = 2;          % whether second-order or first-order system
        degree              % degree of (polynomial) nonlinearity of the rhs of the dynamical system
        nKappa              % Fourier Series expansion order for Fext
        kappas =   []       % matrix with all kappas in its rows
        
        spectrum = []       % data structure constructed by linear_spectral_analysis method
        Options = DSOptions()

    end
    
    methods
        %% SET methods
        function set.A(obj,A)
            obj.A = A;
            set(obj,'order',1); % since second-order system is assumed by default
        end
        
        function set.fnl(obj,fnl)
            % sets nonlinearity in second order form in multi-index format
            if iscell(fnl) %tensor input
                obj.fnl = fnl_tensor2multi(fnl);
            else
                obj.fnl = fnl;
            end
        end

        function set.fnl_semi(obj,fnl_semi)
           % Nonlinear internal force function
                obj.fnl_semi = fnl_semi;
        end
        
        function set.F_semi(obj,F_semi)
            F_semi = cell(numel(F_semi),1);
            for j = 2:numel(F_semi)
                % ensure double output
                if ~isempty(F_semi{j})
                    F_semi{j} = @(invecs) double(F_semi{j}(invecs));
                end
            end
            obj.F_semi = F_semi;
        end
        
        function set.dF_semi(obj,dF_semi)
            dF_semi = cell(numel(dF_semi),1);
            for j = 2:numel(dF_semi)
                % ensure double output
                if ~isempty(dF_semi{j})
                    dF_semi{j} = @(invecs) double(dF_semi{j}(invecs));
                end
            end            
            obj.dF_semi = dF_semi;
        end
        
        function set.fnl_non(obj,fnl_non)
           % Nonlinear internal force function
                obj.fnl_non = fnl_non;
        end        
        
        function set.dfnl_non(obj,dfnl_non)
            % Nonlinear internal force function
            obj.dfnl_non = dfnl_non;
        end
        
        %% GET methods
              
        function A     = get.A(obj)
            
            if obj.order ==1
                
                if ~isempty(obj.K) %Underlying system is second order
                    A = [-obj.K,         sparse(obj.n,obj.n);
                            sparse(obj.n,obj.n),   obj.M];
                else                
                    A = obj.A;
                end
            elseif obj.order == 2
                A = [-obj.K,         sparse(obj.n,obj.n);
                    sparse(obj.n,obj.n),   obj.M];
            end
            
        end
        
        function B     = get.B(obj)
            if obj.order ==1
                
                if isempty(obj.B) && isempty(obj.K)
                    B = speye(obj.N,obj.N);
                elseif ~isempty(obj.K) %Underlying system is second order
                    B = [obj.C,    obj.M;
                        obj.M,  sparse(obj.n,obj.n)];
                else
                    B = obj.B;
                end
                
            elseif obj.order == 2
                
                B = [obj.C,    obj.M;
                    obj.M,  sparse(obj.n,obj.n)];
            end
        end
        
        function BinvA = get.BinvA(obj)
            BinvA = [sparse(obj.n,obj.n), speye(obj.n,obj.n)
                        -obj.M\obj.K,   -obj.M\obj.C];
        end
        
        function n     = get.n(obj)
            n = length(obj.M);
        end
        
        function N     = get.N(obj)
            N = length(obj.A);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GET methods - Internal Forces
        
        % Intrusive nonlinearity
        function F = get.F(obj)
                
            switch obj.Options.notation
                
                case 'tensor'
                    d = length(obj.fnl) + 1;
                    F = cell(1,d);
                    F{1} = sptensor(obj.A);
                    
                    for j = 2:d
                        sizej = obj.N*ones(1,j+1);
                        if isempty(obj.fnl(j-1)) || isempty(obj.fnl(j-1).coeffs)
                            F{j} = sptensor(sizej);
                        else
                            [fnl_t] = multi_index_to_tensor(obj.fnl(j-1).coeffs,obj.fnl(j-1).ind);
                            subsj = fnl_t.subs;
                            valsj = -fnl_t.vals;
                            if obj.order==1
                                valsj = -valsj;
                            end
                            F{j} = sptensor(subsj,valsj,sizej);
                        end
                    end
                    
                case 'multiindex'
                    if obj.order == 1 && ~isempty(obj.F)
                        F = obj.F;
                    else
                        d = length(obj.fnl) + 1;
                        F = repmat(struct('coeffs',[],'ind',[]),1,d);
                        
                        for j = 2:d
                            if isempty(obj.fnl(j-1))
                                F(j) = [];
                                %% following elseif could be removed if inputs are always strictly either first
                                % or second order - not like in bernoulli beam
                            elseif size(obj.fnl(j-1).coeffs,1) == obj.N %fnl already 1st order form, i.e BernoulliBeam
                                F(j).coeffs = obj.fnl(j-1).coeffs;
                                F(j).ind    = obj.fnl(j-1).ind;
                                
                            else % conversion to 1st order form
                                F(j).coeffs = [-obj.fnl(j-1).coeffs;...
                                    sparse(obj.n, size(obj.fnl(j-1).coeffs,2)) ];
                                if obj.n == size(obj.fnl(j-1).ind,2) % No nonlinear damping
                                    F(j).ind = [obj.fnl(j-1).ind.';...
                                        sparse(obj.n, size(obj.fnl(j-1).ind,1)) ].';
                                else %Nonlinear damping
                                    F(j).ind = obj.fnl(j-1).ind;
                                end
                                
                            end
                        end
                    end
                otherwise
                    error('The option should be tensor or multiindex.');
            end
        end
        
        % Semi-intrusive nonlinearity
        function fnl_semi  = get.fnl_semi(obj)
            % sets nonlinearity in second order form in multi-index format
            fnl_semi = obj.fnl_semi;
        end
        
        function dfnl_semi = get.dfnl_semi(obj)
            % sets nonlinearity in second order form in multi-index format
            dfnl_semi = obj.dfnl_semi;
        end
        
        function F_semi    = get.F_semi(obj)
            if obj.order ==1 && ~isempty(obj.F_semi)
                F_semi = obj.F_semi;
                
            elseif obj.order == 2 || (isempty(obj.F_semi) && ~isempty(obj.fnl_semi)) % second order system, or second order and computation carried out first order
                F_semi = set_F_semi(obj);
                
            end
        end
        
        function dF_semi   = get.dF_semi(obj)
            if obj.order ==1 && ~isempty(obj.dF_semi)
                dF_semi = obj.dF_semi;
                
            elseif obj.order == 2 || (isempty(obj.dF_semi) && ~isempty(obj.dfnl_semi)) % second order system, or second order and computation carried out first order
                dF_semi = set_dF_semi(obj);
                
            end
        end
        
        % Non-intrusive nonlinearity
        function fnl_non  = get.fnl_non(obj)
            % sets nonlinearity in second order form in multi-index format
            fnl_non = obj.fnl_non;
        end        

        function dfnl_non = get.dfnl_non(obj)
            % sets nonlinearity in second order form in multi-index format
            dfnl_non = obj.dfnl_non;
        end     
        
        function F_non    = get.F_non(obj)
            if obj.order ==1  && ~isempty(obj.F_non)
                F_non = obj.F_non;
                
            elseif obj.order == 2 || (isempty(obj.F_non) && ~isempty(obj.fnl_non)) % second order system, or second order and computation carried out first order
                F_non = set_F_non(obj);
                
            end
        end

        function dF_non   = get.dF_non(obj)
            if obj.order ==1 && ~isempty(obj.dF_non)
                dF_non = obj.dF_non;
                
            elseif obj.order == 2 || (isempty(obj.dF_non) && ~isempty(obj.dfnl_non)) % second order system, or second order and computation carried out first order
                dF_non = set_dF_non(obj);
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % GET methods - External forcing
        
        function nKappa = get.nKappa(obj)
            nKappa = numel(obj.Fext.data);
        end
        
        function kappas = get.kappas(obj)
            %kappas stored in rows
            sz_kappa = size(obj.Fext.data(1).kappa,2);
            kappas = reshape([obj.Fext.data.kappa],sz_kappa,[]).';
            
        end
        
        function Fext   = get.Fext(obj)
            if obj.order == 1 && ~isempty(obj.Fext)
                Fext = obj.Fext;
            elseif obj.order == 2 || ( isempty(obj.Fext) && ~isempty(obj.fext)) % second order system with second or first order computation
                Fext.data    = set_Fext(obj);
                Fext.epsilon = obj.fext.epsilon;
            end                           
        end
               
        function degree = get.degree(obj)
            degree = 0;
            if ~isempty(obj.A)
                degree = length(obj.F);
            end
        end
        
        %% other methods
        
        [V, D, W] = linear_spectral_analysis(obj)
        
        function add_forcing(obj,f,varargin)
            if ~isfield(f, 'data') % new format
                switch obj.order
                    case 1
                        nn = size(f,1);
                        Kappas = varargin{1};
                        data(1).kappa = Kappas(1);
                        data(2).kappa = Kappas(2);
                        if nn == obj.N                            
                            data(1).F_n_k(1).coeffs = f(:,1);
                            data(1).F_n_k(1).ind    = sparse(1,obj.N);
                            data(2).F_n_k(1).coeffs = f(:,2);
                            data(2).F_n_k(1).ind    = sparse(1,obj.N);
                        elseif nn == obj.N/2 % second order system, first order computation
                            data(1).F_n_k(1).coeffs = [f(:,1);sparse(nn,1)];
                            data(1).F_n_k(1).ind    = sparse(1,obj.N);
                            data(2).F_n_k(1).coeffs = [f(:,2);sparse(nn,1)];
                            data(2).F_n_k(1).ind    = sparse(1,obj.N);
                        else
                            error('Forcing vector has wrong dimension.')
                        end
                        f_ext.data = data;
                        
                    case 2
                        nn = size(f,1);
                        Kappas = varargin{1};
                        data(1).kappa = Kappas(1);
                        data(2).kappa = Kappas(2);
                        data(1).f_n_k(1).coeffs = f(:,1);
                        data(1).f_n_k(1).ind    = sparse(1,nn);
                        data(2).f_n_k(1).coeffs = f(:,2);
                        data(2).f_n_k(1).ind    = sparse(1,nn);
                        f_ext.data = data;
                end
                
                if numel(varargin)>1
                    f_ext.epsilon = varargin{2};
                else
                    f_ext.epsilon = 1;
                end
            else
                f_ext = f;
            end
            
            switch obj.order

                case 1

                    obj.Fext = f_ext;                        

                    if isfield(f_ext,'epsilon')
                        obj.Fext.epsilon = f_ext.epsilon;
                        
                    elseif nargin == 3
                        obj.Fext.epsilon = varargin{1};

                    else
                        obj.Fext.epsilon = 1;

                    end

                case 2                   
                    obj.fext.data = f_ext.data;                        

                    if isfield(f_ext,'epsilon')
                        obj.fext.epsilon = f_ext.epsilon;
                        
                    elseif nargin == 3
                        obj.fext.epsilon = varargin{1};

                    else
                        obj.fext.epsilon = 1;

                    end
            end
        end
        
        fext = compute_fext(obj,t,x)
        Fext = evaluate_Fext(obj,t,z)
        fnl  = compute_fnl(obj,x,xd)
        dfnl = compute_dfnldx(obj,x,xd)
        dfnl = compute_dfnldxd(obj,x,xd)
        Fnl  = evaluate_Fnl(obj,z)
        f    = odefun(obj,t,z)
        [r, drdqdd,drdqd,drdq, c0] = residual(obj, q, qd, qdd, t)
    end
end


function [F_semi]    = set_F_semi(obj)
% Checks input dimensions of vectors and unifies to f_func
F_semi = cell(numel(obj.fnl_semi)+1,1 );

for j = 2:numel(obj.fnl_semi)+1
    %invecs is a N times nl_order dimensional array containing a input
    %vector in each column
    if ~isempty(obj.fnl_semi{j-1})
    F_semi{j} = @(invecs) [-double(obj.fnl_semi{j-1}(invecs)); sparse(obj.n,1)];
    end
end
end

function [dF_semi]   = set_dF_semi(obj)
% Checks input dimensions of vectors and unifies to f_func
dF_semi = cell(numel(obj.dfnl_semi)+1,1 );

for j = 2:numel(obj.dfnl_semi)+1
    %invecs is a N times nl_order dimensional array containing a input
    %vector in each column
    if ~isempty(obj.dfnl_semi{j-1})

            dF_semi{j} = @(invecs) [-double(obj.dfnl_semi{j-1}(invecs)); sparse(obj.n,1)];
    end
end
end

function [F_non]     = set_F_non(obj)

F_non = @(v) [-obj.fnl_non(v); sparse(obj.n,1)];

end

function [dF_non]    = set_dF_non(obj)

dF_non = @(v) [-obj.dfnl_non(v), sparse(obj.n,obj.n); sparse(obj.n,2*obj.n)];

end

function [fnl_multi] = fnl_tensor2multi(fnlTensor)
%Sets second order nonlinear force in multi-index format
d   = length(fnlTensor) + 1;

fnl_multi = repmat(struct('coeffs',[],'ind',[]),1,d-1);

for j = 2:d
    if isempty(fnlTensor{j-1}) || nnz(fnlTensor{j-1}) == 0

    else
        sizej = fnlTensor{j-1}.size;
        subsj = fnlTensor{j-1}.subs;
        valsj = fnlTensor{j-1}.vals;
        tmp = tensor_to_multi_index(sptensor(subsj,valsj,sizej));
        fnl_multi(j-1).coeffs = tmp.coeffs;
        fnl_multi(j-1).ind = tmp.ind;
    end
end
end

function [data]      = set_Fext(obj)
% Creates the data struct from the input second order force
% Structs for storing the coefficients
F_n_k = repmat(struct('coeffs',[],'ind',[]),numel(obj.fext.data(1).f_n_k),1);
data  = repmat(struct('kappa',[],'F_n_k',[]),numel(obj.fext.data),1);

% Fill the structs
for i = 1:numel(obj.fext.data)
    for j = 1:numel(obj.fext.data(i).f_n_k)
        F_n_k(j).coeffs = [obj.fext.data(i).f_n_k(j).coeffs;...
            sparse(obj.n, size(obj.fext.data(i).f_n_k(j).coeffs,2)) ];
        F_n_k(j).ind = [obj.fext.data(i).f_n_k(j).ind.';...
            sparse(obj.n, size(obj.fext.data(i).f_n_k(j).ind,1)) ].';
    end
    data(i).F_n_k = F_n_k;
    data(i).kappa = obj.fext.data(i).kappa;
end
end