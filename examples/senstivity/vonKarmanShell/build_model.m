function [M,C,K,fnl,fext,outdof] = build_model(nDiscretization,E,rho,dampType)

% build unit model
[M0,K0,fnl0,fext,outdof] = build_unit_model(nDiscretization);

% recale matrices and tensors
M = M0*rho;
K = K0*E;
fnl = cell(1,2);
fnl{1} = fnl0{1}*E;
fnl{2} = fnl0{2}*E;

%% Damping matrix
switch dampType
    case 'type1'
        disp('Using Rayleigh damping')
        al = 0.402153037834110;
        be = 8.631167461080214e-06;
        C  = al*M+be*K;

    case 'type2'
        disp('With fixed damping ratio for the first mode');
        lamd = eigs(K,M,1,'smallestabs'); 
        zeta = 0.02;
        om = sqrt(lamd);
        al = zeta*om; 
        be = zeta/om;
        C  = al*M+be*K;
end

end
