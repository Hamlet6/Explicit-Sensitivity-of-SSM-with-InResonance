function [Rk] = Aut_2ndOrder_RedDyn(I,F,THETA,PHI,C,Lambda_K,Lambda, M,Vm,Ym,l,z_k)
% AUT_2NDORDER_REDDYN
%
% This function computes the autonomous reduced dynamics coefficients at
% order k using the second order systems routine. 
%
% [Rk] = AUT_2NDORDER_REDDYN(I,F,THETA,PHI,C,Lambda_K,Lambda, M,Vm,Ym,l,z_k)
%
% F:        multi-index positions that are resonant
% I:        Eigenvalue positions that are resonant
% THETA:    Left Eigenvectors of system
% PHI:      Right Eigenvectors of system
% C:        Damping matrix of system
% Lambda_K:  direct product of multi-indices at order K with master mode
%           vector
% Lambda:   Vector containing master mode eigenvalues
% M:        Mass matrix
% Vm:       Velocity variable RHS of invariance equation
% Ym:       Displacement variable RHS of invariance equation
% l:        dimension of SSM
% z_k:      number of distinct multi-indices at order k
%
% Rk:       order k autonomous reduced dynamics
%
% See also: AUT_2NDORDER_SSM, NONAUT_2NDORDER_REDDYN

Rk = zeros(l,z_k);
% unique multi indices
if any(F)
    [F_un, ~,i_F_un] = unique(F.');
    
    % Loop over multi-indices that lead to resonance
    ii = 1;
    
    for f = F_un
        I_f = I((i_F_un == ii)); % All resonant eigenvalues for this multi - index
        THETA_f = THETA(:,I_f);
        PHI_f   = PHI(:,I_f);
        
        % Set analogous to first order case
        %Rk(I_f,f)= Lambda(I_f) .* -THETA_f'*M*Vm(:,f) + -THETA_f' *Ym(:,f);
        
        % Set analytically to second order case
        Rk(I_f,f) = -( Lambda_K(f).*THETA_f' * M*Vm(:,f) + THETA_f' *(Ym(:,f))) / ( THETA_f' * ( C* PHI_f + M * ((Lambda_K(f) + Lambda(I_f).') .* PHI_f)));
       
        % tbd: Implementation for the case of 1:1 internal resonances!
        
        ii = ii +1;
    end
end
end