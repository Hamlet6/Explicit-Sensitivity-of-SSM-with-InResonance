function pBB = comp_bb(ngrids,E,rho,freqRange)

% create model and setup
[M,K,fnl,~,DM,DK,dfnl2,dfnl3,~,outdof] = build_model_semi(ngrids,E,rho);
al = 0.402153037834110;
be = 8.631167461080214e-06;
DS = DynamicalSystem();
set(DS,'M',M,'C',al*M+be*K,'K',K,'fnl_semi',fnl);
set(DS.Options,'Emax',3,'Nmax',6)
% explict SSM
S = SSM(DS);
resonant_modes = [1 2]; % choose master spectral subspace
S.choose_E(resonant_modes);
[We,Re] = S.explicit_whisker();
% explicit derivative
set(DS,'DM',DM,'DK',DK);
set(DS,'al',@(x) al,'be',@(x) be,'dal',@(x) 0,'dbe',@(x) 0);
set(DS,'Dfnl2_semi',dfnl2,'Dfnl3_semi',dfnl3);
[DW,DR] = S.explicit_senstivity_whisker(We,Re);
pBB = S.BB_sensitivity(We,Re,DW,DR,outdof,freqRange,1,0);

end