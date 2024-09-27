% function [outModel, M, C, K, L, Null, Nullf, ud, outdof, uscale,M_0,C_0,K_0] = build_model_nonIntrusive(alpha)
import com.comsol.model.*
import com.comsol.model.util.*
alpha = 0.1;
model = ModelUtil.create('Model');

model.modelPath('/Users/xuzhenwei/Desktop/Delft/SSMtool-research-NonAndSemiIntr_2_Unify/examples/2DBeam_comsol');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').physics.create('solid', 'SolidMechanics', 'geom1');

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');
model.study('std1').feature('stat').activate('solid', true);

model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('size', [1 0.05]);
model.component('comp1').geom('geom1').run('r1');
model.component('comp1').geom('geom1').run;

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').label('1050 [solid]');
model.component('comp1').material('mat1').info.create('DIN');
model.component('comp1').material('mat1').info('DIN').body('Al99.5');
model.component('comp1').material('mat1').info.create('UNS');
model.component('comp1').material('mat1').info('UNS').body('A91050');
model.component('comp1').material('mat1').info.create('Composition');
model.component('comp1').material('mat1').info('Composition').body('99.5 Al min, 0.05 Cu max, 0.4 Fe max, 0.05 Mg max, 0.05 Mn max, 0.25 Si max, 0.05 V max (wt%)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', 'k_solid_1(T[1/K])[W/(m*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:thermalconductivity', ['Reference: C.Y. Ho, R.W. Powell, P.E. Liley, Journal of Physical and Chemical Reference Data, v1, No. 2, p279 (1972) available online at https://srd.nist.gov/JPCRD/jpcrd7.pdf, same data as for elemental aluminum\nNote: well-annealed with residual resistivity of 0.000594 uohm-cm, error is 2-3% near RT, 3-5% at others, 99.9999%']);
model.component('comp1').material('mat1').propertyGroup('def').set('resistivity', 'res_solid_1(T[1/K])[ohm*m]');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:resistivity', ['Reference: P.D. Desai, H.M. James, C.Y. Ho, Journal of Physical and Chemical Reference Data, v13, No. 4, p1131 (1984) available online at https://srd.nist.gov/JPCRD/jpcrd260.pdf\nNote: data below -233.1C (40.0K) is for Al with a residual resistivity of 1 x 10E-12 ohm-m, not corrected for thermal expansion, purity 99.9% or higher']);
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', '(alpha_solid_1(T[1/K])[1/K]+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha_solid_1(T[1/K])[1/K]-alpha_solid_1(Tempref[1/K])[1/K])/(T-Tempref),d(alpha_solid_1(T[1/K])[1/K],T)))/(1+alpha_solid_1(Tempref[1/K])[1/K]*(Tempref-293[K]))');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:thermalexpansioncoefficient', ['Reference: F.C. Nix, D. MacNair, Physical Review, v60, p597 (1941) and R. Feder, A.S. Norwick, Physical Review, v109, p1959 (1958) and D.F. Gibbons, Physical Review, v112, p136 (1958)\nNote: the reference temperature is 20C (293K), 8% error, same data as elemental Al\nReference temperature: 293.00[K]']);
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'C_solid_1(T[1/K])[J/(kg*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:heatcapacity', ['Reference: B.J. McBride, S. Gordon, M.A. Reno, NASA Technical Paper 3287 (1993)\nNote: same data as elemental aluminum']);
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', 'sigma_solid_1(T[1/K])[S/m]');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:electricconductivity', ['Reference: P.D. Desai, H.M. James, C.Y. Ho, Journal of Physical and Chemical Reference Data, v13, No. 4, p1131 (1984) available online at https://srd.nist.gov/JPCRD/jpcrd260.pdf\nNote: data below -233.1C (40.0K) is for Al with a residual resistivity of 1 x 10E-12 ohm-m, not corrected for thermal expansion, purity 99.9% or higher']);
model.component('comp1').material('mat1').propertyGroup('def').set('HC', 'HC_solid_1(T[1/K])[J/(mol*K)]');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:HC', ['Reference: B.J. McBride, S. Gordon, M.A. Reno, NASA Technical Paper 3287 (1993)\nNote: same data as elemental aluminum']);
model.component('comp1').material('mat1').propertyGroup('def').set('VP', 'VP_solid_1(T[1/K])[Pa]');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:VP', ['Reference: C.B. Alcock, V.P. Itkin and M.K. Horrigan, Canadian Metallurgical Quarterly, v23, p309 (1984)\nNote: same data as elemental Al, 5% error or less']);
model.component('comp1').material('mat1').propertyGroup('def').set('emissivity', 'epsilon(T[1/K])');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:emissivity', ['Reference: K.G. Ramanathan, S.H. Yen and E.A. Estalote, Applied Optics, v11, p2810 (1970)\nNote: electropolished']);
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho_solid_1(T[1/K])[kg/m^3]');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:density', ['Reference: F.C. Nix, D. MacNair, Physical Review, v60, p597 (1941) and R. Feder, A.S. Norwick, Physical Review, v109, p1959 (1958) and D.F. Gibbons, Physical Review, v112, p136 (1958)\nNote: same data as elemental Al']);
model.component('comp1').material('mat1').propertyGroup('def').set('TD', 'TD(T[1/K])[m^2/s]');
model.component('comp1').material('mat1').propertyGroup('def').set('INFO_PREFIX:TD', 'Reference: calculated from the thermal conductivity, density and specific heat');
model.component('comp1').material('mat1').propertyGroup('def').func.create('k_solid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('k_solid_1').set('funcname', 'k_solid_1');
model.component('comp1').material('mat1').propertyGroup('def').func('k_solid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k_solid_1').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('k_solid_1').set('pieces', {'293.0' '933.0' '39.64599+1.684012*T^1-0.005413421*T^2+8.431302E-6*T^3-6.537049E-9*T^4+2.002031E-12*T^5'});
model.component('comp1').material('mat1').propertyGroup('def').func('k_solid_1').label('Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('k_solid_1').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('k_solid_1').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').func.create('res_solid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('res_solid_1').set('funcname', 'res_solid_1');
model.component('comp1').material('mat1').propertyGroup('def').func('res_solid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('res_solid_1').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('res_solid_1').set('pieces', {'1.0' '20.0' '1.091612E-12-1.10726E-13*T^1+3.696901E-14*T^2-2.781934E-15*T^3+1.008733E-16*T^4';  ...
'20.0' '50.0' '-3.313487E-11+7.29041E-12*T^1-4.771551E-13*T^2+1.071535E-14*T^3';  ...
'50.0' '200.0' '1.0965563E-10-3.988929E-11*T^1+1.061978E-12*T^2-2.337666E-15*T^3';  ...
'200.0' '933.0' '-1.037048E-8+1.451201E-10*T^1-8.192563E-14*T^2+6.619834E-17*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('res_solid_1').label('Piecewise 1');
model.component('comp1').material('mat1').propertyGroup('def').func('res_solid_1').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('res_solid_1').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').func.create('alpha_solid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('alpha_solid_1').set('funcname', 'alpha_solid_1');
model.component('comp1').material('mat1').propertyGroup('def').func('alpha_solid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('alpha_solid_1').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('alpha_solid_1').set('pieces', {'20.0' '220.0' '1.371347E-5+7.808536E-8*T^1-2.568882E-10*T^2+3.615726E-13*T^3'; '220.0' '610.0' '5.760185E-6+1.707141E-7*T^1-6.548135E-10*T^2+1.220625E-12*T^3-1.064883E-15*T^4+3.535918E-19*T^5'; '610.0' '933.0' '1.9495E-5+9.630182E-9*T^1+9.462013E-13*T^2'});
model.component('comp1').material('mat1').propertyGroup('def').func('alpha_solid_1').label('Piecewise 2');
model.component('comp1').material('mat1').propertyGroup('def').func('alpha_solid_1').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('alpha_solid_1').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').func.create('C_solid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('C_solid_1').set('funcname', 'C_solid_1');
model.component('comp1').material('mat1').propertyGroup('def').func('C_solid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('C_solid_1').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('C_solid_1').set('pieces', {'100.0' '320.0' '-290.416126+11.1810036*T^1-0.0412540099*T^2+7.11275398E-5*T^3-4.60821994E-8*T^4'; '320.0' '933.0' '595.658507+1.51302896*T^1-0.00207006538*T^2+1.30360846E-6*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('C_solid_1').label('Piecewise 3');
model.component('comp1').material('mat1').propertyGroup('def').func('C_solid_1').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('C_solid_1').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').func.create('sigma_solid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('sigma_solid_1').set('funcname', 'sigma_solid_1');
model.component('comp1').material('mat1').propertyGroup('def').func('sigma_solid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('sigma_solid_1').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('sigma_solid_1').set('pieces', {'1.0' '20.0' '1/(1.008733E-16*T^4-2.781934E-15*T^3+3.696901E-14*T^2-1.107260E-13*T+1.091612E-12)';  ...
'20.0' '50.0' '1/(1.071535E-14*T^3-4.771551E-13*T^2+7.290410E-12*T-3.313487E-11)';  ...
'50.0' '200.0' '1/(-2.337666E-15*T^3+1.061978E-12*T^2-3.988929E-11*T+1.0965563E-10)';  ...
'200.0' '933.0' '1/(6.619834E-17*T^3-8.192563E-14*T^2+1.451201E-10*T-1.037048E-08)'});
model.component('comp1').material('mat1').propertyGroup('def').func('sigma_solid_1').label('Piecewise 4');
model.component('comp1').material('mat1').propertyGroup('def').func('sigma_solid_1').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('sigma_solid_1').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').func.create('HC_solid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('HC_solid_1').set('funcname', 'HC_solid_1');
model.component('comp1').material('mat1').propertyGroup('def').func('HC_solid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('HC_solid_1').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('HC_solid_1').set('pieces', {'100.0' '320.0' '-7.83586214+0.301680207*T^1-0.00111309504*T^2+1.91912758E-6*T^3-1.24336723E-9*T^4'; '320.0' '933.0' '16.0717565+0.0408237901*T^1-5.58534712E-5*T^2+3.51733064E-8*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('HC_solid_1').label('Piecewise 5');
model.component('comp1').material('mat1').propertyGroup('def').func('HC_solid_1').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('HC_solid_1').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').func.create('VP_solid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('VP_solid_1').set('funcname', 'VP_solid_1');
model.component('comp1').material('mat1').propertyGroup('def').func('VP_solid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('VP_solid_1').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('VP_solid_1').set('pieces', {'293.0' '933.0' '(exp((-1.73420000e+04/T-7.92700000e-01*log10(T)+1.23398100e+01)*log(10.0)))*1.33320000e+02'});
model.component('comp1').material('mat1').propertyGroup('def').func('VP_solid_1').label('Piecewise 6');
model.component('comp1').material('mat1').propertyGroup('def').func('VP_solid_1').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('VP_solid_1').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').func.create('epsilon', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('epsilon').set('funcname', 'epsilon');
model.component('comp1').material('mat1').propertyGroup('def').func('epsilon').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('epsilon').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('epsilon').set('pieces', {'175.0' '700.0' '-0.006159012+1.209023E-4*T^1-1.728543E-7*T^2+1.274369E-10*T^3'});
model.component('comp1').material('mat1').propertyGroup('def').func('epsilon').label('Piecewise 7');
model.component('comp1').material('mat1').propertyGroup('def').func('epsilon').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('epsilon').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').func.create('rho_solid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('rho_solid_1').set('funcname', 'rho_solid_1');
model.component('comp1').material('mat1').propertyGroup('def').func('rho_solid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho_solid_1').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('rho_solid_1').set('pieces', {'20.0' '130.0' '2734.317-0.02751647*T^1+0.001016054*T^2-1.700864E-5*T^3+5.734155E-8*T^4'; '130.0' '933.0' '2736.893-0.006011681*T^1-7.012444E-4*T^2+1.3582E-6*T^3-1.367828E-9*T^4+5.177991E-13*T^5'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho_solid_1').label('Piecewise 8');
model.component('comp1').material('mat1').propertyGroup('def').func('rho_solid_1').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('rho_solid_1').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').func.create('TD', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func('TD').set('funcname', 'TD');
model.component('comp1').material('mat1').propertyGroup('def').func('TD').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('TD').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('def').func('TD').set('pieces', {'293.0' '933.0' '5.850406E-5+4.372409E-7*T^1-1.684332E-9*T^2+2.848785E-12*T^3-2.339518E-15*T^4+7.401871E-19*T^5'});
model.component('comp1').material('mat1').propertyGroup('def').func('TD').label('Piecewise 9');
model.component('comp1').material('mat1').propertyGroup('def').func('TD').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('def').func('TD').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('def').addInput('strainreferencetemperature');
model.component('comp1').material('mat1').propertyGroup.create('ThermalExpansion', 'Thermal expansion');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').set('dL', '(dL_solid_1(T[1/K])-dL_solid_1(Tempref[1/K]))/(1+dL_solid_1(Tempref[1/K]))');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').set('INFO_PREFIX:dL', ['Reference: F.C. Nix, D. MacNair, Physical Review, v60, p597 (1941) and R. Feder, A.S. Norwick, Physical Review, v109, p1959 (1958) and D.F. Gibbons, Physical Review, v112, p136 (1958)\nNote: the reference temperature is 20C (293K), same data as elemental Al\nReference temperature: 293.00[K]']);
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').set('alphatan', 'CTE(T[1/K])[1/K]');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').set('INFO_PREFIX:alphatan', ['Reference: F.C. Nix, D. MacNair, Physical Review, v60, p597 (1941) and D.F. Gibbons, Physical Review, v112, p136 (1958)\nNote: same data as elemental Al']);
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func.create('dL_solid_1', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('dL_solid_1').set('funcname', 'dL_solid_1');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('dL_solid_1').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('dL_solid_1').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('dL_solid_1').set('pieces', {'20.0' '188.0' '-0.004116601-4.000347E-6*T^1+5.370388E-8*T^2+3.714324E-10*T^3-1.45073E-12*T^4'; '188.0' '933.0' '-0.00631208932+2.156284E-5*T^1-4.744254E-9*T^2+1.811015E-11*T^3-7.336673E-15*T^4'});
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('dL_solid_1').label('Piecewise');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('dL_solid_1').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('dL_solid_1').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func.create('CTE', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('CTE').set('funcname', 'CTE');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('CTE').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('CTE').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('CTE').set('pieces', {'20.0' '79.0' '-3.317274E-6+3.068688E-7*T^1-1.004816E-8*T^2+1.724768E-10*T^3-8.846061E-13*T^4'; '79.0' '230.0' '-2.288239E-5+6.674915E-7*T^1-4.402622E-9*T^2+1.455358E-11*T^3-1.910622E-14*T^4'; '230.0' '900.0' '1.243109E-5+5.050772E-8*T^1-5.806556E-11*T^2+3.014305E-14*T^3'});
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('CTE').label('Piecewise 1');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('CTE').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').func('CTE').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('ThermalExpansion').addInput('strainreferencetemperature');
model.component('comp1').material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.component('comp1').material('mat1').propertyGroup('Enu').set('youngsmodulus', 'E(T[1/K])[Pa]');
model.component('comp1').material('mat1').propertyGroup('Enu').set('INFO_PREFIX:youngsmodulus', ['Reference: room temperature value from ASM Handbook, v2, 10th edition, ASM International (1992) and temperature dependence from M. Lalpoor, D.G. Eskin, L. Katgerman, Metallurgical and Materials Transactions A, v40A, No. 13, p3304 (2009) and S.C. Sharma, Metallurgical and Materials Transactions A, v31A, No. 3, p773 (2000) and K. Sakai, A. Matsumuro, M. Senoo, Journal of Materials Science, v31, No. 12, p3309 (1996) and E.R. Naimon, H.M. Ledbetter, W.F. Weston, Journal of Materials Science, v10, p1309 (1975) and R.B. McLellan, T. Ishikawa, Journal of Physics and Chemistry of Solids, v48, No. 7, p603 (1987) and H.J. Stokes, Scientific Instruments, v37, p117 (1960)\nNote: uncertainty at -273.1C (0.0K) is 5%, at 500C (773K) it is 10%']);
model.component('comp1').material('mat1').propertyGroup('Enu').set('poissonsratio', 'nu(T[1/K])');
model.component('comp1').material('mat1').propertyGroup('Enu').set('INFO_PREFIX:poissonsratio', ['Reference: temperature dependence from M. Lalpoor, D.G. Eskin, L. Katgerman, Metallurgical and Materials Transactions A, v40A, No. 13, p3304 (2009) and S.C. Sharma, Metallurgical and Materials Transactions A, v31A, No. 3, p773 (2000) and K. Sakai, A. Matsumuro, M. Senoo, Journal of Materials Science, v31, No. 12, p3309 (1996) and E.R. Naimon, H.M. Ledbetter, W.F. Weston, Journal of Materials Science, v10, p1309 (1975) and R.B. McLellan, T. Ishikawa, Journal of Physics and Chemistry of Solids, v48, No. 7, p603 (1987) and H.J. Stokes, Scientific Instruments, v37, p117 (1960)\nNote: uncertainty at -273.1C (0.0K) is 5%, at 500C (773K) it is 10%']);
model.component('comp1').material('mat1').propertyGroup('Enu').func.create('E', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('Enu').func('E').set('funcname', 'E');
model.component('comp1').material('mat1').propertyGroup('Enu').func('E').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('Enu').func('E').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('Enu').func('E').set('pieces', {'0.0' '773.0' '7.770329E10+2036488.0*T^1-189160.7*T^2+425.2931*T^3-0.3545736*T^4'});
model.component('comp1').material('mat1').propertyGroup('Enu').func('E').label('Piecewise');
model.component('comp1').material('mat1').propertyGroup('Enu').func('E').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('Enu').func('E').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('Enu').func.create('nu', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('Enu').func('nu').set('funcname', 'nu');
model.component('comp1').material('mat1').propertyGroup('Enu').func('nu').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('Enu').func('nu').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('Enu').func('nu').set('pieces', {'0.0' '773.0' '0.3238668+3.754548E-6*T^1+2.213647E-7*T^2-6.565023E-10*T^3+4.21277E-13*T^4+3.170505E-16*T^5'});
model.component('comp1').material('mat1').propertyGroup('Enu').func('nu').label('Piecewise 1');
model.component('comp1').material('mat1').propertyGroup('Enu').func('nu').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('Enu').func('nu').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('Enu').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup.create('KG', 'Bulk modulus and shear modulus');
model.component('comp1').material('mat1').propertyGroup('KG').set('G', 'mu(T[1/K])[Pa]');
model.component('comp1').material('mat1').propertyGroup('KG').set('INFO_PREFIX:G', ['Reference: temperature dependence from M. Lalpoor, D.G. Eskin, L. Katgerman, Metallurgical and Materials Transactions A, v40A, No. 13, p3304 (2009) and S.C. Sharma, Metallurgical and Materials Transactions A, v31A, No. 3, p773 (2000) and K. Sakai, A. Matsumuro, M. Senoo, Journal of Materials Science, v31, No. 12, p3309 (1996) and E.R. Naimon, H.M. Ledbetter, W.F. Weston, Journal of Materials Science, v10, p1309 (1975) and R.B. McLellan, T. Ishikawa, Journal of Physics and Chemistry of Solids, v48, No. 7, p603 (1987) and H.J. Stokes, Scientific Instruments, v37, p117 (1960)\nNote: uncertainty at -273.1C (0.0K) is 5%, at 500C (773K) it is 10%']);
model.component('comp1').material('mat1').propertyGroup('KG').set('K', 'kappa(T[1/K])[Pa]');
model.component('comp1').material('mat1').propertyGroup('KG').set('INFO_PREFIX:K', ['Reference: temperature dependence from M. Lalpoor, D.G. Eskin, L. Katgerman, Metallurgical and Materials Transactions A, v40A, No. 13, p3304 (2009) and S.C. Sharma, Metallurgical and Materials Transactions A, v31A, No. 3, p773 (2000) and K. Sakai, A. Matsumuro, M. Senoo, Journal of Materials Science, v31, No. 12, p3309 (1996) and E.R. Naimon, H.M. Ledbetter, W.F. Weston, Journal of Materials Science, v10, p1309 (1975) and R.B. McLellan, T. Ishikawa, Journal of Physics and Chemistry of Solids, v48, No. 7, p603 (1987) and H.J. Stokes, Scientific Instruments, v37, p117 (1960)\nNote: errors may be large']);
model.component('comp1').material('mat1').propertyGroup('KG').func.create('mu', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('KG').func('mu').set('funcname', 'mu');
model.component('comp1').material('mat1').propertyGroup('KG').func('mu').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('KG').func('mu').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('KG').func('mu').set('pieces', {'0.0' '773.0' '2.936605E10-121772.0*T^1-70037.91*T^2+160.9718*T^3-0.1368524*T^4'});
model.component('comp1').material('mat1').propertyGroup('KG').func('mu').label('Piecewise');
model.component('comp1').material('mat1').propertyGroup('KG').func('mu').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('KG').func('mu').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('KG').func.create('kappa', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('KG').func('kappa').set('funcname', 'kappa');
model.component('comp1').material('mat1').propertyGroup('KG').func('kappa').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('KG').func('kappa').set('extrap', 'constant');
model.component('comp1').material('mat1').propertyGroup('KG').func('kappa').set('pieces', {'0.0' '773.0' '7.35752E10+1900026.0*T^1-72136.97*T^2+53.47482*T^3'});
model.component('comp1').material('mat1').propertyGroup('KG').func('kappa').label('Piecewise 1');
model.component('comp1').material('mat1').propertyGroup('KG').func('kappa').set('fununit', '');
model.component('comp1').material('mat1').propertyGroup('KG').func('kappa').set('argunit', '');
model.component('comp1').material('mat1').propertyGroup('KG').addInput('temperature');
model.component('comp1').material('mat1').set('family', 'aluminum');

model.study.create('std2');
model.study('std2').create('time', 'Transient');
model.study('std2').feature('time').activate('solid', true);

model.component('comp1').physics('solid').feature('lemm1').create('dmp1', 'Damping', 2);
model.component('comp1').physics('solid').feature('lemm1').feature('dmp1').set('alpha_dM', alpha);
model.component('comp1').physics('solid').feature('lemm1').feature('dmp1').set('beta_dK', 0);
model.component('comp1').physics('solid').create('fix1', 'Fixed', 1);
model.component('comp1').physics('solid').feature('fix1').selection.set([1]);
model.component('comp1').physics('solid').create('pl1', 'PointLoad', 0);
model.component('comp1').physics('solid').feature('pl1').selection.set([4]);

model.component('comp1').mesh('mesh1').run;
model.component('comp1').mesh('mesh1').autoMeshSize(3);
model.component('comp1').mesh('mesh1').run;

model.study('std1').feature('stat').set('geometricNonlinearity', true);
model.study('std2').feature('time').set('geometricNonlinearity', true);

model.sol.create('sol1');
model.sol('sol1').study('std1');

model.study('std1').feature('stat').set('notlistsolnum', 1);
model.study('std1').feature('stat').set('notsolnum', '1');
model.study('std1').feature('stat').set('listsolnum', 1);
model.study('std1').feature('stat').set('solnum', '1');

model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').feature('st1').set('study', 'std1');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').feature('aDef').set('cachepattern', true);
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature('fc1').set('termonres', 'auto');
model.sol('sol1').feature('s1').feature('fc1').set('reserrfact', 1000);
model.sol('sol1').feature('s1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol1').feature('s1').feature('fc1').set('termonres', 'auto');
model.sol('sol1').feature('s1').feature('fc1').set('reserrfact', 1000);
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.sol('sol1').attach('std1');

model.sol('sol1').runAll;


model.sol.create('sol2');
model.sol('sol2').study('std2');

model.study('std2').feature('time').set('notlistsolnum', 1);
model.study('std2').feature('time').set('notsolnum', '1');
model.study('std2').feature('time').set('listsolnum', 1);
model.study('std2').feature('time').set('solnum', '1');

model.sol('sol2').create('st1', 'StudyStep');
model.sol('sol2').feature('st1').set('study', 'std2');
model.sol('sol2').feature('st1').set('studystep', 'time');
model.sol('sol2').create('v1', 'Variables');
model.sol('sol2').feature('v1').feature('comp1_u').set('scalemethod', 'manual');
model.sol('sol2').feature('v1').feature('comp1_u').set('scaleval', '1e-2*1.0012492197250393');
model.sol('sol2').feature('v1').set('control', 'time');
model.sol('sol2').create('t1', 'Time');
model.sol('sol2').feature('t1').set('tlist', 'range(0,0.1,1)');
model.sol('sol2').feature('t1').set('plot', 'off');
model.sol('sol2').feature('t1').set('plotfreq', 'tout');
model.sol('sol2').feature('t1').set('probesel', 'all');
model.sol('sol2').feature('t1').set('probes', {});
model.sol('sol2').feature('t1').set('probefreq', 'tsteps');
model.sol('sol2').feature('t1').set('atolglobalvaluemethod', 'factor');
model.sol('sol2').feature('t1').set('endtimeinterpolation', true);
model.sol('sol2').feature('t1').set('timemethod', 'genalpha');
model.sol('sol2').feature('t1').set('control', 'time');
model.sol('sol2').feature('t1').feature('aDef').set('cachepattern', true);
model.sol('sol2').feature('t1').create('fc1', 'FullyCoupled');
model.sol('sol2').feature('t1').feature('fc1').set('linsolver', 'dDef');
model.sol('sol2').feature('t1').feature.remove('fcDef');
model.sol('sol2').attach('std2');
model.sol('sol2').runAll;
%%

%% get MCK from time-dependent problem
MA = mphmatrix(model ,'sol2', ...
 'Out', {'K','D','E','Null','Nullf','ud','uscale'});
M_0 = MA.E; C_0 = alpha*MA.E; K_0 = MA.K; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud; uscale = MA.uscale;

M = Nullf'*M_0*Null;
C = Nullf'*C_0*Null;
% K = Nullf'*K_0*Null;

%%

model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 300000 0]);

ML = mphmatrix(model ,'sol1','Out', {'L','K','ud','Nullf'});
L = ML.Nullf'*(ML.L-ML.K*ML.ud);
outdof = find(L);

model.component('comp1').physics('solid').feature('pl1').set('Fp', [0 0 0]);
%%

MA = mphmatrix(model ,'sol1', ...
 'Out', {'K','Null','Nullf','ud','uscale'});
K_0 = MA.K; Null = MA.Null;
Nullf = MA.Nullf; ud = MA.ud; uscale = MA.uscale;
K = Nullf'*K_0*Null;


outModel = model;

% end
