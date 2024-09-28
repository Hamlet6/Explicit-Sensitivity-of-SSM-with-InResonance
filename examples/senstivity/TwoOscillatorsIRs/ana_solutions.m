function [W,R,DW,DR] = ana_solutions(m,b1,zeta,varargin)
% this the encoding of the derivation results shown in
% analytic_derivation.mlx. 

om    = 1/sqrt(m);
phi   = [1/sqrt(m);0];
lamd  = (-zeta+1i*sqrt(1-zeta^2))/sqrt(m);
rlamd = real(lamd);

%% linear terms
% solutions
W10 = [phi; lamd*phi];
% derivatives
dlamddm  = -lamd/(2*m);
dphidm   = [-1/(2*m^1.5);0];
dw102dm  = 2*lamd*[-1/(2*m^1.5);0];
dW10  = [dphidm zeros(2); dw102dm zeros(2)];
dlamd = [dlamddm; 0; 0];
drlamddm = real(dlamddm);
dlamdtdm = 2*dlamddm+conj(dlamddm);
%% quadratic terms
tmp1 = 8*lamd^2*m+6*zeta*lamd*sqrt(m)+2+(8*m^2*lamd+4*zeta*m^1.5)*dlamddm;
tmp2 = 4*lamd^2*m+4*zeta*lamd*sqrt(m)+2;
w201 = [0; -1/m/tmp2];
W20  = [w201; 2*lamd*w201];
dw201dm  = [0; tmp1/tmp2^2/m^2];
dw202dm  = 2*dlamddm*w201+2*lamd*dw201dm;
dW20 = [dw201dm zeros(2); dw202dm zeros(2)];
tmp3 = 4*rlamd^2*m+4*zeta*rlamd*sqrt(m)+2;
w111 = [0; -2/m/tmp3];
W11  = [w111; 2*real(lamd)*w111];
tmp4 = 8*rlamd^2*m+6*zeta*rlamd*sqrt(m)+2+(8*m^2*rlamd+4*zeta*m^1.5)*real(dlamddm);
dw111dm = [0; 2*tmp4/m^2/tmp3^2];
dw112dm = 2*real(dlamddm)*w111+2*rlamd*dw111dm;
dW11 = [dw111dm zeros(2); dw112dm zeros(2)];

%% cubic terms
tmp5 = 9*lamd^2*m+6*zeta*lamd*sqrt(m)+1;
w301 = -[1/m^1.5/tmp5*(1-b1/tmp2); 0];
W30  = [w301; 3*lamd*w301];
tmp6 = 22.5*lamd^2*m^1.5+12*lamd*zeta*m+1.5*sqrt(m)+(18*lamd*m^2.5+6*zeta*m^2)*dlamddm;
g30 = tmp6/m^3/tmp5^2;
h30 = b1*(tmp6*tmp2+m^1.5*tmp5*(4*lamd^2+2*zeta*lamd/sqrt(m)+(8*lamd*m+4*zeta*sqrt(m))*dlamddm))/(m^3*tmp5^2*tmp2^2);
dw301dm = [g30-h30; 0];
dw302dm = 3*dlamddm*w301+3*lamd*dw301dm;
dw301db1 = [1/m^1.5/tmp5/tmp2; 0];
dw302db1 = 3*lamd*dw301db1;
dW30 = [dw301dm dw301db1 zeros(2,1); dw302dm dw302db1 zeros(2,1)];
kapa = -1i/(2*om*sqrt(1-zeta^2));
r21  = -kapa/m^2*(3-b1/tmp2-2*b1/tmp3);
dkapadm = 1i/(2*om^2*sqrt(1-zeta^2))*(-0.5/m^1.5);
tmp7 = 12*lamd^2*m^2+10*lamd*zeta*m^1.5+4*m+(8*lamd*m^3+4*zeta*m^2.5)*dlamddm;
tmp8 = 12*rlamd^2*m^2+10*rlamd*zeta*m^1.5+4*m+(8*rlamd*m^3+4*zeta*m^2.5)*real(dlamddm);
g21  = tmp7/m^4/tmp2^2;
h21  = tmp8/m^4/tmp3^2;
dr21dm  = dkapadm*r21/kapa-kapa*(-6/m^3+b1*g21+2*b1*h21);
dr21db1 = -kapa/m^2*(-1/tmp2-2/tmp3);
dr21 = [dr21dm; dr21db1; 0];
tmp9 = 1/m^1.5*(-1/tmp2-2/tmp3);
tmp0 = dr21db1/(2*rlamd);
lamdt = 2*lamd+conj(lamd);
tmpt = lamdt^2*m+2*lamdt*zeta*sqrt(m)+1;
f211 = 1/m^1.5*(3-b1/tmp2-2*b1/tmp3); varrho = r21/(2*rlamd);
w211 = [-f211/tmpt-varrho/sqrt(m); 0];
w212 = [-lamdt*f211/tmpt-lamd*varrho/sqrt(m); 0];
W21  = [w211; w212];
df211dm     = -1.5/m^2.5*(3-b1/tmp2-2*b1/tmp3)+...
    1/m^1.5*(b1*(4*lamd^2+2*lamd*zeta/sqrt(m))/tmp2^2)+...
    1/m^1.5*(2*b1*(4*rlamd^2+2*rlamd*zeta/sqrt(m))/tmp3^2);
df211dlamd  = 1/m^1.5*b1*(8*lamd*m+4*zeta*sqrt(m))/tmp2^2; 
df211drlamd = 2/m^1.5*b1*(8*rlamd*m+4*zeta*sqrt(m))/tmp3^2; 
a21 = -(df211dm+df211dlamd*dlamddm+df211drlamd*drlamddm)*tmpt;
b21 = -(lamdt^2+lamdt*zeta/sqrt(m)+(2*lamdt*m+2*zeta*sqrt(m))*dlamdtdm)*f211;
c21 = lamdt*a21-dlamdtdm*f211*tmpt;
d21 = lamdt*b21;
dvarrhodm = dr21dm/(2*rlamd)-r21/(2*rlamd^2)*drlamddm;
dw211dm = (a21-b21)/tmpt^2-(dvarrhodm*sqrt(m)-0.5*varrho/sqrt(m))/m;
dw211dm = [dw211dm; 0];
dw212dm = (c21-d21)/tmpt^2-((dvarrhodm*lamd+dlamddm*varrho)*sqrt(m)-0.5*lamd*varrho/sqrt(m))/m;
dw212dm = [dw212dm; 0];
dw211db1 = [-tmp9/tmpt-tmp0/sqrt(m); 0];
dw212db1 = [-lamdt*tmp9/tmpt-lamd*tmp0/sqrt(m); 0];
dW21 = [dw211dm dw211db1 zeros(2,1); dw212dm dw212db1 zeros(2,1)];

%% record output
W  = struct();  R  = struct();
DW = struct();  DR = struct();
W.W10   = W10;
W.W20   = W20;
W.W11   = W11;
W.W30   = W30;
W.W21   = W21;
R.lamd  = lamd;
R.gamma = r21;
DW.DW10 = dW10;
DW.DW20 = dW20;
DW.DW11 = dW11;
DW.DW30 = dW30;
DW.DW21 = dW21;
DR.Dlamd = dlamd;
DR.Dgamma = dr21;

% dW10,dW20,dW11,dW30,dW21,dlamd,dr21
%% non-autonomous part
if numel(varargin)>0
    b2   = varargin{1};
    Om   = varargin{2};
    TMP1 = 2*(Om^2*m-2*zeta*sqrt(m)*Om*1i-1);
    TMP2 = 2*(Om^2*m-2*zeta*sqrt(m)*Om*1i-2);
    a01 = kapa*(2*zeta*sqrt(m)+lamd*m+1i*m*Om)-m;
    b01 = Om^2*m-2*zeta*sqrt(m)*Om*1i-1;
    c01 = 1i*Om*(kapa*(2*zeta*sqrt(m)+lamd*m)-m)-m*kapa*Om^2;
    x01 = a01/TMP1;
    x02 = -b1*b2/TMP2;
    x03 = kapa/2+c01/TMP1;
    x04 = -b1*b2*Om*1i/TMP2;
    x0  = [x01;x02;x03;x04];
    W.x0 = x0;
    R.ftilde = 0.5*sqrt(m)*kapa;
    % derivative of x0
    da01dm = dkapadm*(2*zeta*sqrt(m)+lamd*m+1i*m*Om)+kapa*(zeta/sqrt(m)+dlamddm*m+lamd+1i*Om)-1;
    db01dm = (Om^2-zeta*Om*1i/sqrt(m));
    dx01dm = [(da01dm*b01-a01*db01dm)/(2*b01^2); b1*b2*db01dm/(TMP2^2/2)];
    dx01db1 = [0; -b2/TMP2];
    dx01db2 = [0; -b1/TMP2];
    dc01dm = 1i*Om*(dkapadm*(2*zeta*sqrt(m)+lamd*m)+...
        kapa*(zeta/sqrt(m)+dlamddm*m+lamd)-1)-kapa*Om^2-m*dkapadm*Om^2;
    dx02dm  = [0.5*dkapadm+(dc01dm*b01-c01*db01dm)/(2*b01^2); b1*b2*db01dm*Om*1i/(TMP2^2/2)];
    dx02db1 = [0; -b2*Om*1i/TMP2];
    dx02db2 = [0; -b1*Om*1i/TMP2];
    dx0 = [dx01dm dx01db1 dx01db2;
        dx02dm dx02db1 dx02db2];
    DW.Dx0 = dx0;
    DR.Dftilde = [kapa/(4*sqrt(m))+sqrt(m)*dkapadm/2; 0; 0];
end
