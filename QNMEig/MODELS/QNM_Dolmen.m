function out = model
%
% dolmen_PRA2020_multipole.m
%
% Model exported on Jan 15 2020, 10:59 by COMSOL 5.4.0.225.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('E:\Post_doc1\Post_doc1\Programs\QNMEig_toolbox_multipole\COMSOL model');

model.label('QNM_Dolmen.mph');

model.comments(['QNMEig sphere\n\n']);

model.param.set('omegap_Ag', '1.3659e16 [rad/s]', 'Material: Plasma frequency of Lorentz-Drude Permittivity for Silver');
model.param.set('gamma_Ag', '0.0023*omegap_Ag', 'Material: Damping frequency of Lorentz-Drude Permittivity for Silver');
model.param.set('omega0_Ag', '0 [rad/s]', 'Material: Resonance frequency of Lorentz-Drude permittivity for Silver');
model.param.set('epsiloninf_Ag', '1', 'Material: Silver permittivity at infinite large frequency');
model.param.set('epsilonb', '1', 'Material: Background permittivity');
model.param.set('rb', '300[nm]', 'Geom: Outermost background radius');
model.param.set('t_pml', '300 [nm]', 'Geom: Air domain size');
model.param.set('r_pml', 'rb+t_pml', 'Geom: The total domain size');
model.param.set('width_doublet', '30[nm]', 'Geom: Width of doublet');
model.param.set('length_doublet', '100[nm]', 'Geom: Length of doublet');
model.param.set('width_gros', '128[nm]', 'Geom: Width of the top rod');
model.param.set('length_gros', '50[nm]', 'Geom: Length of the top rod');
model.param.set('epaisseur', '20[nm]', 'Geom: Thickness');
model.param.set('gap_gros', '30[nm]', 'Geom: The separation between the top rod and the doublet');
model.param.set('gap_doublet', '30[nm]', 'Geom: Gap of the doublet');
model.param.set('lambda_N', '100 [nm]', 'Normalization parameter');
model.param.set('lambda_pml', '600 [nm]', 'Typical working wavelength for which PMLs should work properly');
model.param.set('fap', 'c_const/600[nm]', 'The guessed eigen freqiency');

model.component.create('mod1', false);

model.component('mod1').geom.create('geom1', 3);

model.component('mod1').mesh.create('mesh1');

model.component('mod1').geom('geom1').geomRep('comsol');
model.component('mod1').geom('geom1').create('blk1', 'Block');
model.component('mod1').geom('geom1').feature('blk1').set('pos', {'0' 'length_gros/2+gap_gros+length_doublet/2-(length_gros+gap_gros)/2' '0'});
model.component('mod1').geom('geom1').feature('blk1').set('base', 'center');
model.component('mod1').geom('geom1').feature('blk1').set('size', {'width_gros' 'length_gros' 'epaisseur'});
model.component('mod1').geom('geom1').create('blk2', 'Block');
model.component('mod1').geom('geom1').feature('blk2').set('pos', {'(gap_doublet+width_doublet)/2' '-(length_gros+gap_gros)/2' '0'});
model.component('mod1').geom('geom1').feature('blk2').set('base', 'center');
model.component('mod1').geom('geom1').feature('blk2').set('size', {'width_doublet' 'length_doublet' 'epaisseur'});
model.component('mod1').geom('geom1').create('blk3', 'Block');
model.component('mod1').geom('geom1').feature('blk3').set('pos', {'-(gap_doublet+width_doublet)/2' '-(length_gros+gap_gros)/2' '0'});
model.component('mod1').geom('geom1').feature('blk3').set('base', 'center');
model.component('mod1').geom('geom1').feature('blk3').set('size', {'width_doublet' 'length_doublet' 'epaisseur'});
model.component('mod1').geom('geom1').create('sph1', 'Sphere');
model.component('mod1').geom('geom1').feature('sph1').set('layername', {'Layer 1'});
model.component('mod1').geom('geom1').feature('sph1').setIndex('layer', 't_pml', 0);
model.component('mod1').geom('geom1').feature('sph1').set('r', 'r_pml');
model.component('mod1').geom('geom1').run;

model.component('mod1').selection.create('sel1', 'Explicit');
model.component('mod1').selection('sel1').set([1 2 3 4 8 9 10 11]);
model.component('mod1').selection.create('sel2', 'Explicit');
model.component('mod1').selection('sel2').set([1 2 3 4 5 8 9 10 11]);
model.component('mod1').selection.create('sel3', 'Explicit');
model.component('mod1').selection('sel3').set([6 7 12]);
model.component('mod1').selection('sel1').label('PML');
model.component('mod1').selection('sel2').label('Air background and its attached PML');
model.component('mod1').selection('sel3').label('Silver dolmen');

model.variable.create('var1');
model.variable('var1').set('QNM_omega', '(lambda/(-j)) [rad/s]', 'QNM eigen-angular frequencies');
model.component('mod1').variable.create('var2');
model.component('mod1').variable('var2').set('DP1x', 'epsilon0_const*(emw.epsilonrxx*P1x+emw.epsilonrxy*P1y+emw.epsilonrxz*P1z)');
model.component('mod1').variable('var2').set('DP1y', 'epsilon0_const*(emw.epsilonryx*P1x+emw.epsilonryy*P1y+emw.epsilonryz*P1z)');
model.component('mod1').variable('var2').set('DP1z', 'epsilon0_const*(emw.epsilonrzx*P1x+emw.epsilonrzy*P1y+emw.epsilonrzz*P1z)');
model.component('mod1').variable('var2').selection.named('sel3');
model.component('mod1').variable.create('var3');
model.component('mod1').variable('var3').set('QN', '2*intAll((emw.Ex*emw.Dx+emw.Ey*emw.Dy+emw.Ez*emw.Dz)*pml1.detInvT)+intMetal((emw.Ex*emw.Dx+emw.Ey*emw.Dy+emw.Ez*emw.Dz)*pml1.detInvT)*fdisp', 'QNM normalization');
model.component('mod1').variable('var3').set('fdisp', '-2*omegap_Ag^2/(QNM_omega^2-omega0_Ag^2-i*QNM_omega*gamma_Ag)+QNM_omega*omegap_Ag^2*(2*QNM_omega-i*gamma_Ag)/(QNM_omega^2-omega0_Ag^2-i*QNM_omega*gamma_Ag)^2');
model.component('mod1').variable.create('var12');
model.component('mod1').variable('var12').set('omegap1', 'omegap_Ag');
model.component('mod1').variable('var12').set('gamma1', 'gamma_Ag');
model.component('mod1').variable('var12').set('omega01', 'omega0_Ag');
model.component('mod1').variable('var12').set('epsiloninf', 'epsiloninf_Ag');
model.component('mod1').variable('var12').selection.named('sel3');

model.component('mod1').view('view1').hideObjects.create('hide1');
model.component('mod1').view('view1').hideObjects.create('hide2');

model.component('mod1').material.create('mat1', 'Common');
model.component('mod1').material.create('mat2', 'Common');
model.component('mod1').material('mat1').selection.named('sel3');
model.component('mod1').material('mat2').selection.named('sel2');

model.component('mod1').cpl.create('intop1', 'Integration');
model.component('mod1').cpl.create('intop2', 'Integration');
model.component('mod1').cpl('intop1').selection.all;
model.component('mod1').cpl('intop2').selection.set([6 7 12]);

model.component('mod1').coordSystem.create('pml1', 'PML');
model.component('mod1').coordSystem('pml1').selection.named('sel1');

model.component('mod1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.component('mod1').physics('emw').create('weak1', 'WeakContribution', 3);
model.component('mod1').physics('emw').feature('weak1').selection.set([6 7 12]);
model.component('mod1').physics.create('w1', 'WeakFormPDE', 'geom1');
model.component('mod1').physics('w1').field('dimensionless').field('P1');
model.component('mod1').physics('w1').field('dimensionless').component({'P1x' 'P1y' 'P1z'});
model.component('mod1').physics('w1').prop('Units').set('DependentVariableQuantity', 'electricfield');
model.component('mod1').physics('w1').selection.named('sel3');

model.component('mod1').mesh('mesh1').create('ftet1', 'FreeTet');
model.component('mod1').mesh('mesh1').create('ftet3', 'FreeTet');
model.component('mod1').mesh('mesh1').create('swe1', 'Sweep');
model.component('mod1').mesh('mesh1').feature('ftet1').selection.named('sel3');
model.component('mod1').mesh('mesh1').feature('ftet1').create('size1', 'Size');
model.component('mod1').mesh('mesh1').feature('ftet3').selection.geom('geom1', 3);
model.component('mod1').mesh('mesh1').feature('ftet3').selection.set([5]);
model.component('mod1').mesh('mesh1').feature('ftet3').create('size1', 'Size');
model.component('mod1').mesh('mesh1').feature('swe1').selection.geom('geom1', 3);
model.component('mod1').mesh('mesh1').feature('swe1').selection.set([1 2 3 4 8 9 10 11]);
model.component('mod1').mesh('mesh1').feature('swe1').create('dis1', 'Distribution');

model.capeopen.label('Thermodynamics Package');

model.component('mod1').variable('var2').label('Displacement fields associated with auxiliary fields 1');
model.component('mod1').variable('var3').label('QNM normalization');
model.component('mod1').variable('var12').label('Lorentz-Drude parameters for silver domain');

model.component('mod1').view('view1').set('renderwireframe', true);
model.component('mod1').view('view1').hideObjects('hide1').init(1);
model.component('mod1').view('view1').hideObjects('hide1').set('sph1(1)', [3 4 7 8 11 12 14 15 18 19 20 21 23 24 25 26]);
model.component('mod1').view('view1').hideObjects('hide2').init(0);
model.component('mod1').view('view1').hideObjects('hide2').set('sph1(1)', [5 6 7 8]);

model.component('mod1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('mod1').material('mat1').propertyGroup('def').set('electricconductivity', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.component('mod1').material('mat1').propertyGroup('def').set('relpermittivity', {'epsiloninf_Ag' '0' '0' '0' 'epsiloninf_Ag' '0' '0' '0' 'epsiloninf_Ag'});
model.component('mod1').material('mat2').propertyGroup('def').set('relpermittivity', {'epsilonb' '0' '0' '0' 'epsilonb' '0' '0' '0' 'epsilonb'});
model.component('mod1').material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('mod1').material('mat2').propertyGroup('def').set('electricconductivity', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});

model.component('mod1').cpl('intop1').label('All-domain integration');
model.component('mod1').cpl('intop1').set('opname', 'intAll');
model.component('mod1').cpl('intop2').label('Metal-domain integration 1');
model.component('mod1').cpl('intop2').set('opname', 'intMetal');

model.component('mod1').coordSystem('pml1').set('ScalingType', 'Spherical');
model.component('mod1').coordSystem('pml1').set('wavelengthSourceType', 'userDefined');
model.component('mod1').coordSystem('pml1').set('typicalWavelength', 'lambda_pml');

model.component('mod1').physics('emw').prop('MeshControl').set('EnableMeshControl', false);
model.component('mod1').physics('emw').prop('MeshControl').set('SizeControlParameter', 'UserDefined');
model.component('mod1').physics('emw').prop('AnalysisMethodology').set('MethodologyOptions', 'Robust');
model.component('mod1').physics('emw').feature('weak1').set('weakExpression', 'mu0_const*QNM_omega^2*(test(emw.Ex)*DP1x+test(emw.Ey)*DP1y+test(emw.Ez)*DP1z)*pml1.detInvT');
model.component('mod1').physics('emw').feature('weak1').label('Auxililary Field Weak Contribution');
model.component('mod1').physics('w1').prop('ShapeProperty').set('shapeFunctionType', 'shcurl');
model.component('mod1').physics('w1').feature('wfeq1').set('weak', {'1/lambda_N^2*((test(P1x)*P1x+test(P1y)*P1y+test(P1z)*P1z)*(QNM_omega^2-j*gamma_Ag*QNM_omega-omega0_Ag^2)/omegap_Ag^2+(test(P1x)*emw.Ex+test(P1y)*emw.Ey+test(P1z)*emw.Ez))'; '0'; '0'});

model.component('mod1').mesh('mesh1').feature('ftet1').label('scatterer');
model.component('mod1').mesh('mesh1').feature('ftet1').feature('size1').set('custom', 'on');
model.component('mod1').mesh('mesh1').feature('ftet1').feature('size1').set('hmax', '40[nm]/6');
model.component('mod1').mesh('mesh1').feature('ftet1').feature('size1').set('hmaxactive', true);
model.component('mod1').mesh('mesh1').feature('ftet1').feature('size1').set('hmin', '40[nm]/12');
model.component('mod1').mesh('mesh1').feature('ftet1').feature('size1').set('hminactive', true);
model.component('mod1').mesh('mesh1').feature('ftet3').label('air');
model.component('mod1').mesh('mesh1').feature('ftet3').feature('size1').set('custom', 'on');
model.component('mod1').mesh('mesh1').feature('ftet3').feature('size1').set('hmax', '0.6E-7');
model.component('mod1').mesh('mesh1').feature('ftet3').feature('size1').set('hmaxactive', true);
model.component('mod1').mesh('mesh1').feature('swe1').set('facemethod', 'quadlegacy52a');
model.component('mod1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('eig', 'Eigenfrequency');
model.study('std1').feature('eig').set('activate', {'emw' 'on' 'w1' 'on'});

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol('sol1').feature('e1').create('i1', 'Iterative');
model.sol('sol1').feature('e1').feature('i1').create('mg1', 'Multigrid');
model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('pr').create('sv1', 'SORVector');
model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('po').create('sv1', 'SORVector');

model.result.create('pg1', 'PlotGroup3D');
model.result.create('pg2', 'PlotGroup3D');
model.result('pg1').create('mslc1', 'Multislice');
model.result('pg1').create('arwv1', 'ArrowVolume');
model.result('pg2').create('slc1', 'Slice');

model.study('std1').feature('eig').set('neigs', 2);
model.study('std1').feature('eig').set('neigsactive', true);
model.study('std1').feature('eig').set('eigunit', 'Hz');
model.study('std1').feature('eig').set('shift', 'fap');
model.study('std1').feature('eig').set('discretization', {'emw' 'physics' 'w1' 'physics'});

model.sol('sol1').attach('std1');
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('neigs', 2);
model.sol('sol1').feature('e1').set('shift', 'fap');
model.sol('sol1').feature('e1').feature('dDef').active(true);
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('e1').feature('i1').label('Suggested Iterative Solver (emw)');
model.sol('sol1').feature('e1').feature('i1').set('linsolver', 'bicgstab');
model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('pr').feature('sv1').set('sorvecdof', {'mod1_E'});
model.sol('sol1').feature('e1').feature('i1').feature('mg1').feature('po').feature('sv1').set('sorvecdof', {'mod1_E'});
model.sol('sol1').runAll;

model.result('pg1').label('Electric Field (emw)');
model.result('pg1').set('looplevel', [2]);
model.result('pg1').set('titletype', 'none');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').set('legendcolor', 'white');
model.result('pg1').feature('mslc1').label('Multislice');
model.result('pg1').feature('mslc1').set('expr', 'emw.normE/abs(sqrt(QN))');
model.result('pg1').feature('mslc1').set('unit', '');
model.result('pg1').feature('mslc1').set('descr', 'emw.normE/abs(sqrt(QN))');
model.result('pg1').feature('mslc1').set('xnumber', '0');
model.result('pg1').feature('mslc1').set('ynumber', '0');
model.result('pg1').feature('mslc1').set('rangecoloractive', true);
model.result('pg1').feature('mslc1').set('rangecolormax', 2.75E16);
model.result('pg1').feature('mslc1').set('colortable', 'Thermal');
model.result('pg1').feature('mslc1').set('resolution', 'fine');
model.result('pg1').feature('mslc1').set('smooth', 'internal');
model.result('pg1').feature('mslc1').set('resolution', 'fine');
model.result('pg1').feature('arwv1').active(false);
model.result('pg1').feature('arwv1').set('expr', {'emw.Ex/emw.normE' 'emw.Ey/emw.normE' 'emw.Ez/emw.normE'});
model.result('pg1').feature('arwv1').set('descr', '');
model.result('pg1').feature('arwv1').set('arrowxmethod', 'coord');
model.result('pg1').feature('arwv1').set('xcoord', 'range(-70[nm],15[nm],70[nm])');
model.result('pg1').feature('arwv1').set('arrowymethod', 'coord');
model.result('pg1').feature('arwv1').set('ycoord', 'range(-90[nm],15[nm],90[nm])');
model.result('pg1').feature('arwv1').set('arrowzmethod', 'coord');
model.result('pg1').feature('arwv1').set('zcoord', 0);
model.result('pg1').feature('arwv1').set('scale', '1e-8');
model.result('pg1').feature('arwv1').set('scaleactive', true);
model.result('pg1').feature('arwv1').set('color', 'white');
model.result('pg2').feature('slc1').set('expr', 'P1x');
model.result('pg2').feature('slc1').set('descr', 'Dependent variable P1, x component');
model.result('pg2').feature('slc1').set('smooth', 'internal');
model.result('pg2').feature('slc1').set('resolution', 'normal');

out = model;
