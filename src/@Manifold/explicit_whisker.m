function [W,R] = explicit_whisker(obj)
% Invariant manifold in the autonomous system limit

% setup
W  = struct();
R  = struct();
Ap = obj.System.A;
Bp = obj.System.B;
Q  = obj.System.F_semi{2};
C  = obj.System.F_semi{3};
midx = obj.System.N;
if ~obj.System.nl_damp && obj.System.order==2 
    midx = midx/2;  % 2nd order fct handle only takes displacement as input
end
ispropDamp = obj.System.Options.RayleighDamping; % damping type
Lambda_E = obj.E.spectrum;
assert(numel(Lambda_E)==2,'current implementation only supports two-dimensional SSMs');
idx = 1;
if imag(Lambda_E(1))<0; idx=2; end
Lambda_E = Lambda_E(idx);
vE  = obj.E.basis(:,idx); vEbar = conj(vE);
uE  = obj.E.adjointBasis(:,idx);

switch obj.System.order
    case 2
        %% second order computation
        fprintf('\n second order explict solution \n')
        % setup
        phi = vE(1:obj.System.n); psi = conj(uE(1:obj.System.n));
        M   = obj.System.M;
        C   = obj.System.C;
        K   = obj.System.K;
        f2  = obj.System.fnl_semi{1};
        f2  = @(x1,x2) f2({x1;x2});
        f3  = obj.System.fnl_semi{2};
        f3  = @(x1,x2,x3) f3({x1;x2;x3});       
        if ispropDamp % proportional damping
            assert(isreal(phi),'eigenvector is not real');
            tmp = phi'*M*phi;
            phi = phi/sqrt(tmp);
            om  = sqrt(phi'*K*phi);
            zet = phi'*C*phi/(2*om);
            assert(abs(Lambda_E-(-zet*om+1i*om*sqrt(1-zet^2)))<1e-6*om,'eigenvalue does not match');
            kapa = -1i/(2*om*sqrt(1-zet^2));
            psi = kapa*phi;
        end
        W.psi = psi;
        % linear terms
        W10 = [phi;Lambda_E*phi];    
        solver = obj.Options.solver;

        % quadratic order terms
        A20  = 4*Lambda_E^2*M+2*Lambda_E*C+K;
        if obj.System.nl_damp
            W201 = -solveinveq(A20,f2(vE,vE),solver);
        else
            W201 = -solveinveq(A20,f2(phi,phi),solver);
        end
        W20  = [W201; 2*Lambda_E*W201];
        A11  = 4*(real(Lambda_E))^2*M+2*real(Lambda_E)*C+K;
        if obj.System.nl_damp
            W111 = -solveinveq(A11,f2(vE,vEbar)+f2(vEbar,vE),solver);
        else
            W111 = -2*solveinveq(A11,f2(phi,phi),solver);
        end
        W11  = [W111; 2*real(Lambda_E)*W111];
        
        % cubic terms
        A30  = 9*Lambda_E^2*M+3*Lambda_E*C+K;
        if obj.System.nl_damp
            f30  = f2(vE,W20)+f2(W20,vE)+f3(vE,vE,vE);
        else
            f30  = f2(phi,W201)+f2(W201,phi)+f3(phi,phi,phi);
        end
        W301 = -solveinveq(A30,f30,solver);
        W30  = [W301; 3*Lambda_E*W301];
        lamdt = 2*Lambda_E+conj(Lambda_E);
        A21  = lamdt^2*M+lamdt*C+K;
        if obj.System.nl_damp
            f21  = f2(vEbar,W20)+f2(W20,vEbar)+f2(vE,W11)+f2(W11,vE)+...
                f3(vEbar,vE,vE)+f3(vE,vEbar,vE)+f3(vE,vE,vEbar);
        else
            f21  = f2(phi,W201)+f2(W201,phi)+f2(phi,W111)+f2(W111,phi)+3*f3(phi,phi,phi);
        end
        r21  = -psi.'*f21;
        tmp1 = r21/(2*real(Lambda_E));
        tmp2 = solveinveq(A21,f21,solver);
        W21  = [-tmp2-tmp1*phi; -lamdt*tmp2-tmp1*Lambda_E*phi];

        % leading-order of non-autonomous part
        fext = obj.System.fext;
        if ~isempty(fext)
            ncol = size(fext,2);
            if ncol>1
                assert(ncol==2,'only a pair of forcing is allowed');
                fext = 2*fext(:,1);
            end
            Om  = obj.System.Omega;
            AOm = [-K-1j*Om*C -1j*Om*M; -1j*Om*M M];
            ftilde = 0.5*psi.'*fext;
            bOm = [(C+Lambda_E*M)*phi*ftilde-0.5*fext; M*phi*ftilde];
            x0  = solveinveq(AOm,bOm,solver);
            W.x0 = x0;
            R.ftilde = ftilde;
        end

    case 1
        %% first order computation
        % Linear terms
        W10 = vE;
        
        % quadratic order terms
        solver = obj.Options.solver;
        W20 = solveinveq(2*Lambda_E*Bp-Ap,Q({vE(1:midx);vE(1:midx)}),solver);
        W11 = solveinveq(2*real(Lambda_E)*Bp-Ap,...
            Q({vE(1:midx);vEbar(1:midx)})+Q({vEbar(1:midx);vE(1:midx)}),solver);
        
        % cubic order terms
        W30 = solveinveq(3*Lambda_E*Bp-Ap,Q({vE(1:midx);W20(1:midx)})+...
            Q({W20(1:midx);vE(1:midx)})+C({vE(1:midx);vE(1:midx);vE(1:midx)}),solver);
        tmp = Q({vEbar(1:midx);W20(1:midx)})+Q({W20(1:midx);vEbar(1:midx)})+...
            Q({vE(1:midx);W11(1:midx)})+Q({W11(1:midx);vE(1:midx)})+...
            C({vEbar(1:midx);vE(1:midx);vE(1:midx)})+...
            C({vE(1:midx);vEbar(1:midx);vE(1:midx)})+...
            C({vE(1:midx);vE(1:midx);vEbar(1:midx)});
        r21 = uE'*tmp;
        W21 = solveinveq((2*Lambda_E+conj(Lambda_E))*Bp-Ap,tmp-Bp*vE*r21,solver);
end

W.W10 = W10;
W.W20 = W20;
W.W11 = W11;
W.W30 = W30;
W.W21 = W21;
R.r21 = r21;
R.Lambda_E = Lambda_E;
end

