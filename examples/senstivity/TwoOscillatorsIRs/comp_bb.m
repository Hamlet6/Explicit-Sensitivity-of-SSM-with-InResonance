function pBB = comp_bb(m,b1,b2,zeta,freqRange,outdof)

% create model and setup
[mass,damp,stiff,fnl,~,DM,DK,dfnl2,dfnl3,~] = build_model_sens_direct(m,b1,b2,zeta);
DS = DynamicalSystem();
set(DS,'M',mass,'C',damp,'K',stiff,'fnl_semi',fnl);
set(DS.Options,'Emax',1,'Nmax',2)
[~,~,~] = DS.linear_spectral_analysis();
% explict SSM
S = SSM(DS);
S.choose_E([1 2]);
[We,Re] = S.explicit_whisker();
% explicit derivative
set(DS,'DM',DM,'DK',DK);
set(DS,'al',@(x) 2*0.01*x,'be',@(x) 0,'dal',@(x) 2*0.01,'dbe',@(x) 0);
set(DS,'Dfnl2_semi',dfnl2,'Dfnl3_semi',dfnl3);
[DW,DR] = S.explicit_senstivity_whisker(We,Re);
pBB = S.BB_sensitivity(We,Re,DW,DR,outdof,freqRange,1,0);

end