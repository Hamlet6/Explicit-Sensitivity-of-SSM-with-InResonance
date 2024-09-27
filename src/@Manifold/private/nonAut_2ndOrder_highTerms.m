function [W1,R1] = nonAut_2ndOrder_highTerms(obj,W0,R0,W1,R1,reduced_data)
% NONAUT_2NDORDER_HIGHTERMS This function computes the high order,
%
% parametrisation-coordinate dependent parts of the non-autonomous SSM using
% the second order system computation routine.
%
% [W1,R1] = NONAUT_2NDORDER_HIGHTERMS(obj,W0,R0,W1,R1,reduced_data)
%
% obj:      SSM class object
% W0:       autonomous SSM coefficients
% R0:       autonomous RD coefficients
% W1:       non-autonomous SSM coefficients
% R1:       non-autonomous RD coefficients
% red_data: data struct containing necessary information for computation
%
% W1:       non-autonomous SSM coefficients
% R1:       non-autonomous RD coefficients
%
% See also: NONAUT_1STORDER_HIGHTERMS, NONAUT_2NDORDER_WHISKER


%% Solving for coefficients with k>0
% The coefficient equation for this case reads
%
% $$\sum_{i=1}^{2n}     \underbrace{     \bigg(       \mathbf{(A)}_{bi}
% -        \mathbf{B}_{bi}         \big[            \sum_{j=1}^l k_j  \lambda_{j}
% +            i\langle \mathbf{\Omega}, \mathbf{\eta} \rangle         \big]
% \bigg)     }_{:= (\mathcal{L}_{\mathbf{k},\mathbf{\eta}})_{bi}}    U^i_{\mathbf{k},\mathbf{\eta}}\\=
% \sum_{i=1}^{2n} \mathbf{B}_{bi}\sum_{j=1}^l    \bigg[          \sum_{\mathbf{m},
% \mathbf{u}\in \mathbb{N}^l , \ \mathbf{m+u} - \mathbf{\hat{e}}_j = \mathbf{k}}
% m_j        T^i_{\mathbf{m}}        Q^j_{\mathbf{u},\mathbf{\eta}}        +
% \sum_{\mathbf{m,u} \in \mathbb{N}^l, \ |\mathbf{m}|<k \ \ \mathbf{m+u} - \hat{\mathbf{e}}_j
% = \mathbf{k}} m_j U^i_{\mathbf{m},\mathbf{\eta}} P^j_{\mathbf{u}}    \bigg]
% \\  \ \ \ -    \sum_{\mathbf{n}\in \mathbb{N}^{2n}, |\mathbf{n}|<k}
% F^b_{\mathbf{n},\mathbf{\eta}} \pi_{\mathbf{n,k}}-    \sum_{\mathbf{n}\in \mathbb{N}
% ^{2n}, \ \            |\mathbf{n}| \geq 2}       G^b_{\mathbf{n}}\sigma_{\mathbf{k},
% \mathbf{n}, \mathbf{\eta}}$$

% Get autonomous coefficients and composition coefficients in rev. lex.
% ordering
order = reduced_data.order;
kappas = reduced_data.kappas;
l      = reduced_data.l;
data                 = reduced_data;

data.Lambda_M_vector = obj.E.spectrum;
data.THETA           = obj.E.adjointBasis(1:(data.N/2),:);
data.PHI             = obj.E.basis(1:(data.N/2),:);

data.Mass            = obj.System.M;
data.Damp            = obj.System.C;
data.Stiff           = obj.System.K;

FEXT                = obj.System.Fext;
[NL,reduced_data]   = computationMode(obj,reduced_data); % Check if computation is intrusive / uses potential
reduced_data.sym    = obj.System.F_semi_sym; %if nonlinear forces are symmetric

if reduced_data.intrusion
    [W0,R0,reduced_data.H] = get_autonomous_coeffs(W0,R0,true);
else
    [W0,R0,~] = get_autonomous_coeffs(W0,R0,false);
end

%% Perform Nonautonomous Calculation
% We loop over all orders of spatial multi-indices. Within that there is a loop
% over all the frequency multi-indices.


for k = 1:order
    soltic = tic;
    
    % Nl and Fext for all harmonics
    nltic = tic;
    [FGs] = nonAut_Fext_plus_Fnl(NL,FEXT,reduced_data,k,W0,W1,size(kappas,2));
    nltime = toc(nltic);
    
    eqtime = 0;
    mixtime = 0;
    
    for i = 1:size(kappas,2)
        %% Calculating the RHS
        
        %Forcing and nonlinearity terms
        FG = FGs(i).val;
        
        % Mixed Terms
        mixtic = tic;
        [WR] = nonAut_W1R0_plus_W0R1(reduced_data,k,W0,W1(i),R0,R1(i));
        mixtime = mixtime + toc(mixtic);
        %% Find resonant terms
        
        [ev_idx, multi_idx,Lambda_K]   = nonAut_resonant_terms(k,kappas(:,i),reduced_data,'k'); % F contains multi-index pos.
        
        % save memory
        if ~any(FG)
            FG = sparse(FG);
        end
        if ~any(WR)
            WR = sparse(WR);
        end
        
        % Data to pass into function
        data.i           = i;
        data.k           = k;
        data.I           = ev_idx;
        data.F           = multi_idx;
        data.Lambda_K    = Lambda_K;
        
        eqtic = tic;
        [R1_ik,W1_ik] = nonAut_2ndOrder_SolveInvEq(FG,WR,data);
        eqtime = eqtime + toc(eqtic);
        
        [R1,W1] = nonAut_assembleCoefficients(W1,W1_ik,R1,R1_ik,i,k,l);

        
    end
    
    soltime = toc(soltic);
    
    % save information about runtimes
    obj.solInfoNonAut.timeEstimate(k+1) = obj.solInfoNonAut.timeEstimate(k+1) + soltime;
    obj.solInfoNonAut.nlTime(k+1)       = obj.solInfoNonAut.nlTime(k+1) + nltime;
    obj.solInfoNonAut.mixTime(k+1)      = obj.solInfoNonAut.mixTime(k+1) + mixtime;
    obj.solInfoNonAut.eqTime(k+1)       = obj.solInfoNonAut.eqTime(k+1) + eqtime;
end

end

function [NL,data] = computationMode(obj,data)

data.intrusion = false;
switch obj.System.Options.Intrusion
    case 'none'
        data.mode = 'NonIntrusiveF';
        
        % check input dimensionality of function handles
        if obj.System.order == 2
            data.nl_damp = obj.System.nl_damp;
        end
        
        data.nl_ord = 3;
        NL = obj.System.dF_non;
    case 'semi'
        data.mode = 'SemiIntrusiveF';
        
        % check input dimensionality of function handles
        if obj.System.order == 2
            data.nl_damp = obj.System.nl_damp;
        end
        
        data.nl_ord = numel(obj.System.dF_semi);
        NL = obj.System.dF_semi;
        
    case 'full'
        
        data.nl_ord    = numel(obj.System.F);
        data.mode = 'IntrusiveF';
        data.intrusion = true;
        
        NL = obj.System.F;
end

end


function [W0,R0,H]               = get_autonomous_coeffs(W0,R0,intrusion)
% Sets up the autonomous coefficients used in nonautonomous computation

%These quantities are all in lexicographic ordering, calculations are carried out in reverse
%lexicographic ordering. This is accounted for below.

W0 = coeffs_lex2revlex(W0,'TaylorCoeff');
R0 = coeffs_lex2revlex(R0,'TaylorCoeff');

%composition coefficients of power series
if intrusion
    [H] = get_composition_coeffs(W0);
else
    H = [];
end

end

function [H]          = get_composition_coeffs(W0)
% This function reconstructs the composition coefficients for the computed
% SSM coefficients
%W_0 input in rev-lexicographic ordering, outputs H in rev-lexicographic ordering
data.ordering = 'revlex';

H = cell(1,numel(W0));
H{1} = W0(1).coeffs;

for k = 2:numel(W0)
    data.k = k;
    H{k} = coeffs_composition(W0,H,data);
end
end
