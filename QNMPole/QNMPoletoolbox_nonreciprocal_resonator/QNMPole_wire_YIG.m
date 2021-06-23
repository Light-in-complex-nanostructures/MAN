function out = model
%
% QNMPole_wire_YIG.m
%
% Model exported on Jun 7 2021, 12:48 by COMSOL 5.4.0.225.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('E:\Post_doc1\Post_doc1\co-workers\Philippe Lalanne\non-reciprocity\data\4 released');

model.label('QNMPole_wire_YIG.mph');

model.param.set('Lair', 'a*1.5', 'Geom: air background');
model.param.set('r', '0.35*a', 'Geom: rod radius');
model.param.set('Lpml', 'a/4', 'Geom: PML length');
model.param.set('a', '26[mm]', 'Geom: Yang bing');
model.param.set('sym_factor', '1', 'Geom: symmetry factor');
model.param.set('epsrinf', '15', 'Material: epsilonrinf');
model.param.set('murinf', '1', 'Material: muinf');
model.param.set('epsilonb', '1', 'Material: epsilonb');
model.param.set('lambda_pml', 'c_const/freqg', 'Material: typical absorbing wavelength of PMLs');
model.param.set('freqg', '(8.8466)[GHz]', 'Freq: guessed eigenfrequency');
model.param.set('xdipole', 'a/2', 'WL: x position of the dipole');
model.param.set('ydipole', 'a/5', 'WL: y position of the dipole');
model.param.set('omega', '11e9*2*pi[rad/s]', 'WL: omega');
model.param.group.create('par2');
model.param('par2').set('omegam', '175[mT]*gamma');
model.param('par2').set('gamma', '2*pi*28[GHz/T]');
model.param('par2').set('omega0', 'gamma*Hs*mu0_const');
model.param('par2').set('Hs', '900[Oe]');
model.param('par2').set('alpha', '3e-4');
model.param.group.create('par3');
model.param('par3').set('w', '2*pi*1.7[GHz]');
model.param('par3').set('du_1', 'omegam*(i*alpha*w+omega0+omegam)');
model.param('par3').set('low', '((1+alpha^2)*w^2-2*i*alpha*w*(omega0+omegam)-(omega0+omegam)^2)');
model.param('par3').set('b', 'w*omegam/low/murinf');
model.param('par3').set('d', '(du_1/low+1)/murinf');
model.param('par3').set('mur', 'murinf+murinf*omegam*(omega0+i*alpha*w)/((omega0+i*alpha*w)^2-w^2)');
model.param('par3').set('kappa', 'omegam*w/((omega0+i*alpha*w)^2-w^2)*murinf');
model.param('par2').label('YIG');
model.param('par3').label('Test');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.func.create('an8', 'Analytic');
model.func.create('an2', 'Analytic');
model.func.create('an3', 'Analytic');
model.func.create('an5', 'Analytic');
model.func.create('an6', 'Analytic');
model.func.create('an7', 'Analytic');
model.func('an8').label('epsilon_YIG');
model.func('an8').set('funcname', 'epsilon_YIG');
model.func('an8').set('expr', 'epsrinf');
model.func('an8').set('argunit', 'rad/s');
model.func('an8').set('fununit', '1');
model.func('an2').label('mu_air');
model.func('an2').set('funcname', 'mu_air');
model.func('an2').set('expr', '1');
model.func('an2').set('argunit', 'rad/s');
model.func('an2').set('fununit', '1');
model.func('an2').set('complex', true);
model.func('an3').label('epsilon_air');
model.func('an3').set('funcname', 'epsilon_air');
model.func('an3').set('expr', '1');
model.func('an3').set('argunit', 'rad/s');
model.func('an3').set('fununit', '1');
model.func('an3').set('complex', true);
model.func('an5').label('mur_YIG');
model.func('an5').set('funcname', 'mur_YIG');
model.func('an5').set('expr', 'murinf+omegam*(omega0+i*alpha*x)/((omega0+i*alpha*x)^2-x^2)*murinf');
model.func('an5').set('argunit', 'rad/s');
model.func('an5').set('fununit', '1');
model.func('an5').set('complex', true);
model.func('an5').set('plotargs', {'x' '2[GHz]*2*pi' '5[GHz]*2*pi'});
model.func('an6').label('kappa_YIG');
model.func('an6').set('funcname', 'kappa_YIG');
model.func('an6').set('expr', 'omegam*x/((omega0+i*alpha*x)^2-x^2)*murinf');
model.func('an6').set('argunit', 'rad/s');
model.func('an6').set('complex', true);
model.func('an6').set('plotargs', {'x' '9[GHz]*2*pi' '11[GHz]*2*pi'});
model.func('an7').label('mur_inf');
model.func('an7').set('funcname', 'mur_inf');
model.func('an7').set('expr', '1');
model.func('an7').set('argunit', 'rad/s');
model.func('an7').set('fununit', '1');
model.func('an7').set('complex', true);
model.func('an7').set('plotargs', {'x' '2[GHz]*2*pi' '5[GHz]*2*pi'});

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('r', 'r');
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('r1').set('layername', {'Layer 1'});
model.component('comp1').geom('geom1').feature('r1').setIndex('layer', 'Lpml', 0);
model.component('comp1').geom('geom1').feature('r1').set('layerleft', true);
model.component('comp1').geom('geom1').feature('r1').set('layerright', true);
model.component('comp1').geom('geom1').feature('r1').set('layertop', true);
model.component('comp1').geom('geom1').feature('r1').set('size', {'Lair+Lpml*2' 'Lair+Lpml*2'});
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').feature('pt1').set('p', {'xdipole' 'ydipole'});
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').selection.create('sel1', 'Explicit');
model.component('comp1').selection('sel1').set([10]);
model.component('comp1').selection.create('sel2', 'Explicit');
model.component('comp1').selection('sel2').set([1 2 3 4 5 6 7 8 9]);
model.component('comp1').selection.create('sel3', 'Explicit');
model.component('comp1').selection('sel3').set([1 2 3 4 6 7 8 9]);
model.component('comp1').selection('sel1').label('sca');
model.component('comp1').selection('sel2').label('background+pml');
model.component('comp1').selection('sel3').label('PML');

model.component('comp1').material.create('mat1', 'Common');

model.component('comp1').cpl.create('intop1', 'Integration');
model.component('comp1').cpl('intop1').selection.named('sel1');

model.component('comp1').coordSystem.create('pml1', 'PML');
model.component('comp1').coordSystem('pml1').selection.named('sel3');

model.component('comp1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.component('comp1').physics('emw').create('wee2', 'WaveEquationElectric', 2);
model.component('comp1').physics('emw').feature('wee2').selection.named('sel2');
model.component('comp1').physics('emw').create('lco1', 'LineCurrentOutOfPlane', 0);
model.component('comp1').physics('emw').feature('lco1').selection.set([13]);

model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri2').selection.set([5 10]);
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').create('size2', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').selection.named('sel1');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size2').selection.set([5]);
model.component('comp1').mesh('mesh1').feature('map1').create('dis1', 'Distribution');
model.component('comp1').mesh('mesh1').feature('map1').feature('dis1').selection.set([6 12 15 18]);

model.component('comp1').view('view1').axis.set('xmin', -0.04350161924958229);
model.component('comp1').view('view1').axis.set('xmax', 0.04350161924958229);
model.component('comp1').view('view1').axis.set('ymin', -0.04035469517111778);
model.component('comp1').view('view1').axis.set('ymax', 0.04035469517111778);

model.component('comp1').cpl('intop1').label('scaV');
model.component('comp1').cpl('intop1').set('opname', 'scaV');

model.component('comp1').coordSystem('pml1').set('wavelengthSourceType', 'userDefined');
model.component('comp1').coordSystem('pml1').set('typicalWavelength', 'lambda_pml');

model.component('comp1').physics('emw').prop('components').set('components', 'outofplane');
model.component('comp1').physics('emw').feature('wee1').set('epsilonr', {'epsilon_YIG(omega)'; '0'; '0'; '0'; 'epsilon_YIG(omega)'; '0'; '0'; '0'; 'epsilon_YIG(omega)'});
model.component('comp1').physics('emw').feature('wee1').set('mur', {'mur_YIG(omega)'; 'kappa_YIG(omega)*1i'; '0'; '-kappa_YIG(omega)*1i'; 'mur_YIG(omega)'; '0'; '0'; '0'; 'mur_inf(omega)'});
model.component('comp1').physics('emw').feature('wee1').label('Wave Equation, Electric sca');
model.component('comp1').physics('emw').feature('wee2').set('epsilonr', {'epsilon_air(omega)'; '0'; '0'; '0'; 'epsilon_air(omega)'; '0'; '0'; '0'; 'epsilon_air(omega)'});
model.component('comp1').physics('emw').feature('wee2').set('mur', {'mu_air(omega)'; '0'; '0'; '0'; 'mu_air(omega)'; '0'; '0'; '0'; 'mu_air(omega)'});
model.component('comp1').physics('emw').feature('wee2').label('Wave Equation, Electric background');
model.component('comp1').physics('emw').feature('lco1').set('Iop', 2);

model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', '0.0097500');
model.component('comp1').mesh('mesh1').feature('size').set('hmin', '5.8500E-6');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmax', 'r/10');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmin', 5.85E-6);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size2').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size2').set('hmax', 'r/2');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size2').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').run;

model.component('comp1').physics('emw').feature('wee1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('sigma_mat', 'userdef');

model.study.create('std1');
model.study('std1').create('freq', 'Frequency');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.numerical.create('int1', 'IntVolume');
model.result.numerical('int1').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').create('surf1', 'Surface');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
model.sol('sol1').feature('v1').set('cname', {'freq'});
model.sol('sol1').feature('v1').set('clist', {'1[GHz]'});
model.sol('sol1').feature('s1').set('stol', 0.01);
model.sol('sol1').feature('s1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'1[GHz]'});
model.sol('sol1').feature('s1').feature('p1').set('punit', {'GHz'});
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
model.sol('sol1').feature('s1').feature('d1').label('Suggested Direct Solver (emw)');
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').runAll;

model.result('pg1').label('Electric Field (emw)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').label('Surface');
model.result('pg1').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('resolution', 'normal');

out = model;
