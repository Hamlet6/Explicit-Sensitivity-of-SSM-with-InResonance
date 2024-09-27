function [frce,varargout] = comp_frc(ngrids,E,rho,Omega,Q,epsilon,freqRange)

% create model and setup
[M,K,fnl,fext,DM,DK,dfnl2,dfnl3,dfext,outdof] = build_model_semi(ngrids,E,rho);
DS = DynamicalSystem();
al = 0.402153037834110;
be = 8.631167461080214e-06;
set(DS,'M',M,'C',al*M+be*K,'K',K,'fnl_semi',fnl,'fext',fext);
set(DS.Options,'Emax',3,'Nmax',6)
set(DS,'dfext',dfext,'Omega',Omega);
% explict SSM
S = SSM(DS);
set(S.Options,'contribNonAuto',false)
set(S.Options,'solver','backslash');
resonant_modes = [1 2]; % choose master spectral subspace
S.choose_E(resonant_modes);
[We,Re] = S.explicit_whisker();
% explicit derivative
set(DS,'DM',DM,'DK',DK);
set(DS,'al',@(x) al,'be',@(x) be,'dal',@(x) 0,'dbe',@(x) 0);
set(DS,'Dfnl2_semi',dfnl2,'Dfnl3_semi',dfnl3);
[DW,DR] = S.explicit_senstivity_whisker(We,Re);
optdof = outdof; 
[frce,~] = S.po_sensitivity(We,Re,DW,DR,optdof,Q,epsilon);
pFRC = S.FRC_sensitivity(We,Re,DW,DR,optdof,Q,epsilon,freqRange,1,0);
varargout{1} = pFRC;

end