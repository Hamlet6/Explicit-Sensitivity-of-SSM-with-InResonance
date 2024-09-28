function [frce,varargout] = comp_frc(m,b1,b2,zeta,epsilon,Om,varargin)

% create model and setup
[mass,damp,stiff,fnl,fext,DM,DK,dfnl2,dfnl3,dfext] = build_model_sens_direct(m,b1,b2,zeta);
DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl_semi',fnl);
set(DS.Options,'Emax',1,'Nmax',2)
[~,~,~] = DS.linear_spectral_analysis();
set(DS,'fext',fext,'Omega',Om,'dfext',dfext);
% explict SSM
S = SSM(DS);
S.choose_E([1 2]);
[We,Re] = S.explicit_whisker();
% explicit derivative
set(DS,'DM',DM,'DK',DK);
set(DS,'al',@(x) 2*0.01*x,'be',@(x) 0,'dal',@(x) 2*0.01,'dbe',@(x) 0);
set(DS,'Dfnl2_semi',dfnl2,'Dfnl3_semi',dfnl3);
[DW,DR] = S.explicit_senstivity_whisker(We,Re);
optdof = [1 2 3 4]; Q = eye(4); 
[frce,~] = S.po_sensitivity(We,Re,DW,DR,optdof,Q,epsilon);
if numel(varargin)>0
    freqRange = varargin{1};
    pFRC = S.FRC_sensitivity(We,Re,DW,DR,optdof,Q,epsilon,freqRange,1,0);
    varargout{1} = pFRC;
end
end