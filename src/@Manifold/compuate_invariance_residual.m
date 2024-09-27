function y = compuate_invariance_residual(obj,W0,R0,p,type,varargin)
% COMPUTE_INVARIANCE_RESIDUAL This function calcualte the residuals of
% points on the computed SSM. It compute the Euclidean norm of residual of
% the invariance equation for a point on the computed SSM.
%
% y = compuate_invariance_residual(obj,W0,R0,p,type,varargin)
%
% W0 - autonomous part of SSM
% R0 - autonomous part of reduced dynamics on SSM
% p  - parameterized coordinates (dimSSM x npts), it can be a single point
% or a collection points, e.g., a time history
% type - auto/nonauto: for auto type, we take invariance equation for
% autonomous SSM. For nonauto, full invariance equation is considered
% varargin: [W1,R1,] nonautonomous part of SSM and reduced dynamics - needed
% when type = nonauto


% BDW*R=A*W+F(W)+Fext(\Omega t); R=R_auto+R_nonauto

npts = size(p,2);
switch type
    case 'auto'
        % from parameterized p to z in full system
        z = reduced_to_full_traj([],p,W0);
        % left side
        DW = 0; R = 0;
        for j = 1:length(W0)
            DW = DW+expand_multiindex_derivative(W0(j),p); 
            R  = R+expand_multiindex(R0(j),p);
        end        
        lhs = spblkdiag(DW)*R(:);
        lhs = reshape(lhs,obj.dimSystem,npts);
        lhs = obj.System.B*real(lhs);
        % rhs
        rhs = obj.System.A*z+obj.System.evaluate_Fnl(z);
        res = lhs-rhs;
        y = sqrt(sum(res.^2));
        
    case 'nonauto'
        error('not implemented yet');
        W1 = varargin{1}; 
        R1 = varargin{2}; used when high-order nonauto is included
        omega = obj.System.Omega; t = linspace(0,2*pi/omega,npts);
        if obj.System.order==2
            epf = obj.System.fext.epsilon;
        else
            epf = obj.System.Fext.epsilon;
        end
        z = reduced_to_full_traj(t,p,W0,W1,epsilon,omega,epf);  
        
    otherwise 
        error('type should be either auto or nonauto');
end

end