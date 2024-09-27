function varargout = harmonic_balance_NLvib(obj,omegaRange,varargin)
% HARMONIC_BALANCE_NLVIB This function performs periodic orbit
% continuation using the harmonic balance method in the forcing frequency
% range omegaRange
% This function essentially uses the NLvib package.
% optional arguments for NLvib:
% 'nHarmonics': number of harmonics of the forcing frequency to be
%               considered in periodic response approximation
% 'nt': number of equally spaced time intervals over which the
%                periodic response is approximated in time domain
% 'outdof': the degree of freedom at which the FRC is plotted
% 'maxSteps': maximum number of continuation steps for NLvib (5000 by default)
% 'outScale': output scaling (only used for plotting)
% 'stepSize': continuation step size in NLvib
% 'qScale'  : qscl parameter of NLvib
% 'nSample' : sampling rate (only used for plotting)
startHB = tic;
assert(obj.order == 2, 'NLvib can only be used for second-order systems')
n = obj.n;
% read inputs
[nHarmonics, outdof, nt, maxSteps, stepSize, outScale, qScale, nSample] = ...
    parse_inputs(varargin{:});

% get multi-index coefficients
M.coeffs = [];
M.ind = [];

for j = 1:numel(obj.fnl)    
    M.coeffs = cat(2,M.coeffs,obj.fnl(j).coeffs);
    M.ind = cat(1,M.ind,obj.fnl(j).ind);    
end
coefficients = M.coeffs;
powers = M.ind;

% forcing shape vector assuming single harmonic cosine forcing
f_ext = obj.fext.epsilon * real(sum([obj.fext.data(1).f_n_k(1).coeffs,obj.fext.data(2).f_n_k(1).coeffs],2));
system = System_with_PolynomialStiffnessNonlinearity(obj.M,obj.C,obj.K,full(powers),coefficients.',f_ext);

% Analysis parameters
omegaStart = omegaRange(1);      % start frequency
omegaEnd = omegaRange(end);      % end frequency
analysis = 'FRF';

% Initial guess (solution of linearized system)
Q1 = (-omegaStart^2*system.M + 1i*omegaStart*system.D + system.K)\system.Fex1;
y0 = zeros((2*nHarmonics+1)*n,1);
y0(n + (1:2*n)) = [real(Q1); -imag(Q1)];

% Solve and continue w.r.t. Om
options = struct('Dscale',[qScale*ones(size(y0));omegaStart],'stepmax',maxSteps, 'dsmin',stepSize/100);
[X,solInfo] = solve_and_continue(y0,...
    @(X) HB_residual(X,system,nHarmonics,nt,analysis),...
    omegaStart,omegaEnd,stepSize,options);

% save and plot solution
[X, Omega] = deal(X(1:end-1,:), X(end,:));

save('HB.mat', 'X', 'Omega', 'nHarmonics', 'solInfo')

if outdof    
    [varargout{1}, varargout{2}] = plot_HB_sol(X, Omega, nHarmonics,nt,outdof,outScale,nSample);
end
totalComputationTime = toc(startHB);
disp(['Total time spent on FRC computation using harmonic balance via NLvib = ' datestr(datenum(0,0,0,0,0,totalComputationTime),'HH:MM:SS')])

end


function [nHarmonics, outdof, nt, maxSteps, stepSize, outScale, ...
    qScale, nSample] = parse_inputs(varargin)
%% parsing inputs
defaultnHarmonics = 10;
defaultoutdof = 0;
defaultnt = 2^7;
defaultmaxSteps = 5000;
defaultstepSize = 1e-2;
defaultoutScale = 1;
defaultqScale = 7e-1;
defaultnSample = 1;

p = inputParser;
addParameter(p,'nHarmonics',defaultnHarmonics, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'maxSteps',defaultmaxSteps, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'nt',defaultnt, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','positive'}) );
addParameter(p,'outdof',defaultoutdof, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','nonnegative'}) );
addParameter(p,'stepSize',defaultstepSize, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','positive'}) );
addParameter(p,'outScale',defaultoutScale, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','positive'}) );
addParameter(p,'qScale',defaultqScale, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','positive'}) );
addParameter(p,'nSample',defaultnSample, @(x)validateattributes(x, ...
    {'numeric'},{'nonempty','integer','positive'}) );
parse(p,varargin{:});

nHarmonics = p.Results.nHarmonics;
outdof = p.Results.outdof;
nt = p.Results.nt;
maxSteps = p.Results.maxSteps;
stepSize = p.Results.stepSize;
outScale = p.Results.outScale;
qScale = p.Results.qScale;
nSample = p.Results.nSample;
end

function [Adof,Omega] = plot_HB_sol(X, Omega, H, N, dof, outscale, ns)
% This function plots and returns the Amplitude, Frequency values
% at a given degree of freedom in the system.
% Input:   X - NLVib solution (obtained from continuation)
%          H - Number of harmonics
%          N - Number of time samples used in inverse FFT
%          dof - output DOF in the range 1 to n (total number of DOFs)


n = (size(X,1))/(2*H+1);
nOmega = length(Omega);

% Vector recovering deflection at dof
T_tip = sparse(1,n);
T_tip(dof) = 1;

Qtip_HB = kron(eye(2*H+1),T_tip)*X;
Adof = nan(nOmega,1);

for k = 1:ns:nOmega
    Qdum = Qtip_HB(:,k);
    Qfft = zeros(H+1,1);
    Qfft(1)=Qdum(1);
    for j = 1:H
        Qfft(j+1) = Qdum(2*j)+1i*Qdum(2*j+1);
    end
    Adof(k) = outscale*(N/2)*max(abs(ifft(Qfft,N,'symmetric')));
end

hold on
plot(Omega,Adof,'kx','MarkerSize',8,'LineWidth',2,'DisplayName','HB-NLvib')
legend('show')
end