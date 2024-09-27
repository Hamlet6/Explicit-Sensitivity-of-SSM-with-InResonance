function varargout = SSM_sensitivity(obj,W,DW,rhosamp,plotdof,pidx,pfracs)
% FRC_SENSITIVITY This function returns the perturbed SSM when
% the pidx-th parameter is perturbed with pfractions. Here we assume pidx
% is a scalar, yet the pfractions can be a vector.
%
% W:  a structure array with the SSM parameterization
% DW: a structure array with the sensitivity of SSM parameterization
% rhosamp: mesh for rho grid
% plotdof: visualization for SSM
% pidx:    index for perturbed parameters
% pfracs:  perturbation fractions (should be much less than one)

% setup
nfrac = numel(pfracs);
ndof  = numel(plotdof);
zmesh = cell(nfrac,1);
assert(ndof==3,'the number of dofs is not three');
for k=1:nfrac
    % compute perturbed coefficients - expansion coefficients
    Wk = struct();
    Wk.W10 = W.W10+DW.DW10(:,pidx)*pfracs(k);
    Wk.W11 = W.W11+DW.DW11(:,pidx)*pfracs(k);
    Wk.W20 = W.W20+DW.DW20(:,pidx)*pfracs(k);
    Wk.W30 = W.W30+DW.DW30(:,pidx)*pfracs(k);
    Wk.W21 = W.W21+DW.DW21(:,pidx)*pfracs(k);
    % FRC in physical coordinates
    zdofs = ssm_mesh(Wk,rhosamp,plotdof);
    zmesh{k} = zdofs;
end
% compute periodic orbit and its amplitude
varargout{1} = zmesh;
end


function zdofs = ssm_mesh(W,rhosamp,plotdof)
% setup
W10 = W.W10(plotdof,:);
W20 = W.W20(plotdof,:);
W11 = W.W11(plotdof,:);
W30 = W.W30(plotdof,:);
W21 = W.W21(plotdof,:);
% generate grids
[RHO,THETA] = meshgrid(rhosamp,0:2*pi/128:2*pi);
zdof1 = zeros(size(RHO));
zdof2 = zeros(size(RHO));
zdof3 = zeros(size(RHO));
% compute nodes
for k=1:129
    rhok = RHO(k,:);
    thk  = THETA(k,:);
    tmp1 = rhok.*exp(1i*thk);
    tmp2 = rhok.^2.*exp(1i*2*thk);
    tmp3 = rhok.^3.*exp(1i*3*thk);
    tmp4 = rhok.^3.*exp(1i*thk);
    zk = 2*real(W10*tmp1+W20*tmp2+W30*tmp3+W21*tmp4)+W11*(rhok.^2);
    zdof1(k,:) = zk(1,:);
    zdof2(k,:) = zk(2,:);
    zdof3(k,:) = zk(3,:);
end
figure; hold on
h = surf(zdof1,zdof2,zdof3,'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 1.0, ...
    'LineStyle', '-', 'EdgeColor', 0.6*[1 1 1], ...
    'LineWidth', 0.5);
view([1,1,1])
grid on
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
% record output
zdofs = struct();
zdofs.zdof1 = zdof1;
zdofs.zdof2 = zdof2;
zdofs.zdof3 = zdof3;
end
