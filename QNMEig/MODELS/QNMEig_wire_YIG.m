function out = model
%
% QNMEig_wire_YIG.m
%

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('E:\Post_doc1\Post_doc1\co-workers\Philippe Lalanne\non-reciprocity\data\4 released');

model.label('QNMEig_wire_YIG.mph');

model.param.set('Lair', 'a*1.5', 'Geom: air background');
model.param.set('r', '0.35*a', 'Geom: wire radius');
model.param.set('Lpml', 'a/4', 'Geom: PML thickness');
model.param.set('a', '26[mm]', 'Geom: lattice constant');
model.param.set('sym_factor', '1', 'Geom: symmetry factor 1: 0 sym. plane 4: 2 sym. planes');
model.param.set('epsrinf', '15', 'Material: permittivity of YIG');
model.param.set('murinf', '1', 'Material: murinf given in section 1.3.4');
model.param.set('epsilonb', '1', 'Material: permittivity of the background medium');
model.param.set('lambda_pml', 'c_const/freqg', 'Material: typical absorbing wavelength of PMLs');
model.param.set('freqg', '(8.8466)[GHz]', 'Freq: frequency to search for QNMs around');
model.param.set('nomega', '(freqg*2*pi)', 'Freq: normalization factor for the weak form');
model.param.group.create('par2');
model.param('par2').set('omegam1', '175[mT]*gamma1', 'Material: omega_m given in section 1.3.4');
model.param('par2').set('omega01', 'gamma1*Hs1*mu0_const', 'Material: omega_0 given in section 1.3.4');
model.param('par2').set('Hs1', '900[Oe]', 'Material: H0 given in section 1.3.4');
model.param('par2').set('gamma1', '2*pi*28[GHz/T]', 'Material: gamma given in section 1.3.4');
model.param('par2').set('alpha1', '3e-4', 'Material: alpha given in section 1.3.4');
model.param('par2').label('YIG materials');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('selresult', true);
model.component('comp1').geom('geom1').feature('c1').set('color', 'custom');
model.component('comp1').geom('geom1').feature('c1').set('customcolor', [0.3176470696926117 0.48627451062202454 0.9843137264251709]);
model.component('comp1').geom('geom1').feature('c1').set('r', 'r');
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').set('base', 'center');
model.component('comp1').geom('geom1').feature('r1').set('layername', {'Layer 1'});
model.component('comp1').geom('geom1').feature('r1').setIndex('layer', 'Lpml', 0);
model.component('comp1').geom('geom1').feature('r1').set('layerleft', true);
model.component('comp1').geom('geom1').feature('r1').set('layerright', true);
model.component('comp1').geom('geom1').feature('r1').set('layertop', true);
model.component('comp1').geom('geom1').feature('r1').set('size', {'Lair+Lpml*2' 'Lair+Lpml*2'});
model.component('comp1').geom('geom1').run;

model.component('comp1').selection.create('sel1', 'Explicit');
model.component('comp1').selection('sel1').set([10]);
model.component('comp1').selection.create('sel2', 'Explicit');
model.component('comp1').selection('sel2').set([1 2 3 4 5 6 7 8 9]);
model.component('comp1').selection.create('sel3', 'Explicit');
model.component('comp1').selection('sel3').set([1 2 3 4 6 7 8 9]);
model.component('comp1').selection('sel1').label('sca');
model.component('comp1').selection('sel2').label('Air background and its attached PML');
model.component('comp1').selection('sel3').label('PML');

model.component('comp1').variable.create('var2');
model.component('comp1').variable('var2').set('QNM_omega', 'emw.iomega/i', 'Complex eigenfrequency of E^L');
model.component('comp1').variable('var2').set('QNM_omega2', 'emw2.iomega/i', 'Complex eigenfrequency of E^R');
model.component('comp1').variable.create('var5');
model.component('comp1').variable('var5').set('M1x', '0[V/m^2]', 'Auxiliary field M1x in Air domain');
model.component('comp1').variable('var5').set('M1y', '0[V/m^2]', 'Auxiliary field M1y in Air domain');
model.component('comp1').variable('var5').set('N1x', '0[V/m^2]', 'Auxiliary field N1x in Air domain');
model.component('comp1').variable('var5').set('N1y', '0[V/m^2]', 'Auxiliary field N1y in Air domain');
model.component('comp1').variable('var5').selection.named('sel2');
model.component('comp1').variable.create('var3');
model.component('comp1').variable('var3').set('fac1', 'i*alpha1*w2+omega01+omegam1');
model.component('comp1').variable('var3').set('fac3', 'i*w2*omegam1');
model.component('comp1').variable('var3').set('fac2', 'i*alpha1*w2+omega01');
model.component('comp1').variable('var3').set('lowduinvu', '(w2^2-fac1^2)*(w2^2-fac2^2)');
model.component('comp1').variable('var3').set('fac4', '(2*omega01+omegam1)*w2*(1+alpha1^2)');
model.component('comp1').variable('var3').set('fac5', 'i*alpha1*((1+alpha1^2)*w2^2-omega01*(omega01+omegam1))');
model.component('comp1').variable('var3').set('upduinvu11', 'omegam1*(fac4+fac5)');
model.component('comp1').variable('var3').set('upduinvu12', 'i*omegam1*((1+alpha1^2)*w2^2+omega01*(omega01+omegam1))');
model.component('comp1').variable('var3').set('dwudw11', '1/w2+upduinvu11/lowduinvu');
model.component('comp1').variable('var3').set('dwudw12', 'upduinvu12/lowduinvu');
model.component('comp1').variable('var3').set('dwudw21', '-upduinvu12/lowduinvu');
model.component('comp1').variable('var3').set('dwudw22', '1/w2+upduinvu11/lowduinvu');
model.component('comp1').variable('var3').set('dwudwH_x', '(dwudw11*cERx+dwudw12*cERy)/(-i)/mu0_const');
model.component('comp1').variable('var3').set('dwudwH_y', '(dwudw21*cERx+dwudw22*cERy)/(-i)/mu0_const');
model.component('comp1').variable('var3').set('dwudwH_z', 'cERz/w2/(-i)/mu0_const');
model.component('comp1').variable('var3').selection.named('sel1');
model.component('comp1').variable.create('var6');
model.component('comp1').variable('var6').set('dwudwH_x', 'HRx');
model.component('comp1').variable('var6').set('dwudwH_y', 'HRy');
model.component('comp1').variable('var6').set('dwudwH_z', 'HRz');
model.component('comp1').variable('var6').selection.named('sel2');
model.component('comp1').variable.create('var4');
model.component('comp1').variable('var4').set('w', 'emw.iomega/i');
model.component('comp1').variable('var4').set('w2', 'emw2.iomega/i');
model.component('comp1').variable('var4').set('HLx', '1/(-i*w*mu0_const)*(M1x+emw.curlEx*invmuinf)');
model.component('comp1').variable('var4').set('HLy', '1/(-i*w*mu0_const)*(M1y+emw.curlEy*invmuinf)');
model.component('comp1').variable('var4').set('HLz', '1/(-i*w*mu0_const)*(emw.curlEz*invmuinf)');
model.component('comp1').variable('var4').set('invmuinf', '(murinf)^(-1)');
model.component('comp1').variable('var4').set('cERx', 'emw2.curlEx');
model.component('comp1').variable('var4').set('cERy', 'emw2.curlEy');
model.component('comp1').variable('var4').set('cERz', 'emw2.curlEz');
model.component('comp1').variable('var4').set('HRx', '1/(-i*w2*mu0_const)*(N1x+emw2.curlEx*invmuinf)');
model.component('comp1').variable('var4').set('HRy', '1/(-i*w2*mu0_const)*(N1y+emw2.curlEy*invmuinf)');
model.component('comp1').variable('var4').set('HRz', '1/(-i*w2*mu0_const)*(emw2.curlEz*invmuinf)');
model.component('comp1').variable.create('var7');
model.component('comp1').variable('var7').set('mur_YIG', 'murinf+murinf*omegam1*(omega01+i*alpha1*w2)/((omega01+i*alpha1*w2)^2-w2^2)');
model.component('comp1').variable('var7').set('kappa_YIG', 'omegam1*w2/((omega01+i*alpha1*w2)^2-w2^2)*murinf');
model.component('comp1').variable('var7').set('dwkappa_YIG', 'd(kappa_YIG*w2,w2)');
model.component('comp1').variable('var7').set('dmur_YIG', 'd(mur_YIG*w2,w2)');
model.component('comp1').variable('var7').set('dmu11', 'dmur_YIG');
model.component('comp1').variable('var7').set('dmu12', 'dwkappa_YIG*i');
model.component('comp1').variable('var7').set('dmu21', '-dwkappa_YIG*i');
model.component('comp1').variable('var7').set('dmu22', 'dmur_YIG');
model.component('comp1').variable('var7').set('duwHx', 'dmu11*HRx+dmu12*HRy');
model.component('comp1').variable('var7').set('duwHy', 'dmu21*HRx+dmu22*HRy');
model.component('comp1').variable('var7').set('duwHz', 'HRz');
model.component('comp1').variable('var7').selection.named('sel1');
model.component('comp1').variable.create('var8');
model.component('comp1').variable('var8').set('duwHx', 'HRx');
model.component('comp1').variable('var8').set('duwHy', 'HRy');
model.component('comp1').variable('var8').set('duwHz', 'HRz');
model.component('comp1').variable('var8').selection.named('sel2');
model.component('comp1').variable.create('var9');
model.component('comp1').variable('var9').set('lowbd', '((1+alpha1^2)*w2^2-2*i*alpha1*w2*(omega01+omegam1)-(omega01+omegam1)^2)');
model.component('comp1').variable('var9').set('du', 'omegam1*(i*alpha1*w2+omega01+omegam1)');
model.component('comp1').variable('var9').set('ttestHx', '(d*emw2.curlEx+i*b*emw2.curlEy)/(-i)/w2/mu0_const');
model.component('comp1').variable('var9').set('b', 'w2*omegam1/lowbd/murinf');
model.component('comp1').variable('var9').set('d', '(du/lowbd+1)/murinf');
model.component('comp1').variable('var9').selection.named('sel1');

model.component('comp1').coordSystem.create('pml1', 'PML');
model.component('comp1').coordSystem('pml1').selection.named('sel3');

model.component('comp1').physics.create('emw', 'ElectromagneticWaves', 'geom1');
model.component('comp1').physics('emw').create('wee2', 'WaveEquationElectric', 2);
model.component('comp1').physics('emw').feature('wee2').selection.named('sel2');
model.component('comp1').physics('emw').create('weak1', 'WeakContribution', 2);
model.component('comp1').physics('emw').feature('weak1').selection.named('sel1');
model.component('comp1').physics.create('w', 'WeakFormPDE', 'geom1');
model.component('comp1').physics('w').field('dimensionless').field('M1');
model.component('comp1').physics('w').field('dimensionless').component({'M1x' 'M1y'});
model.component('comp1').physics('w').prop('Units').set('DependentVariableQuantity', 'none');
model.component('comp1').physics('w').prop('Units').set('CustomDependentVariableUnit', 'V/m^2');
model.component('comp1').physics('w').selection.named('sel1');
model.component('comp1').physics.create('emw2', 'ElectromagneticWaves', 'geom1');
model.component('comp1').physics('emw2').create('wee2', 'WaveEquationElectric', 2);
model.component('comp1').physics('emw2').feature('wee2').selection.named('sel2');
model.component('comp1').physics('emw2').create('weak1', 'WeakContribution', 2);
model.component('comp1').physics('emw2').feature('weak1').selection.named('sel1');
model.component('comp1').physics.create('w2', 'WeakFormPDE', 'geom1');
model.component('comp1').physics('w2').field('dimensionless').field('N1');
model.component('comp1').physics('w2').field('dimensionless').component({'N1x' 'N1y'});
model.component('comp1').physics('w2').prop('Units').set('DependentVariableQuantity', 'none');
model.component('comp1').physics('w2').prop('Units').set('CustomDependentVariableUnit', 'V/m^2');
model.component('comp1').physics('w2').selection.named('sel1');

model.component('comp1').mesh('mesh1').create('ftri2', 'FreeTri');
model.component('comp1').mesh('mesh1').create('ftri3', 'FreeTri');
model.component('comp1').mesh('mesh1').create('map1', 'Map');
model.component('comp1').mesh('mesh1').feature('ftri2').selection.named('sel1');
model.component('comp1').mesh('mesh1').feature('ftri2').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('ftri3').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('ftri3').selection.set([5]);
model.component('comp1').mesh('mesh1').feature('ftri3').create('size1', 'Size');
model.component('comp1').mesh('mesh1').feature('map1').create('dis1', 'Distribution');
model.component('comp1').mesh('mesh1').feature('map1').feature('dis1').selection.set([6 12 15 18]);

model.component('comp1').variable('var2').label('Eigenfrequency');
model.component('comp1').variable('var5').label('M=0; N=0;');
model.component('comp1').variable('var3').label('QN_dispersion');
model.component('comp1').variable('var6').label('QN_nondispersion');
model.component('comp1').variable('var4').label('Magnetic field');
model.component('comp1').variable('var7').label('duwH');
model.component('comp1').variable('var8').label('duwH_background');
model.component('comp1').variable('var9').active(false);
model.component('comp1').variable('var9').label('b and d');

model.component('comp1').view('view1').set('showselection', false);
model.component('comp1').view('view1').set('showmaterial', true);
model.component('comp1').view('view1').axis.set('xmin', -0.03177777677774429);
model.component('comp1').view('view1').axis.set('xmax', 0.03177777677774429);
model.component('comp1').view('view1').axis.set('ymin', -0.029478957876563072);
model.component('comp1').view('view1').axis.set('ymax', 0.029478957876563072);

model.component('comp1').coordSystem('pml1').set('wavelengthSourceType', 'userDefined');
model.component('comp1').coordSystem('pml1').set('typicalWavelength', 'lambda_pml');

model.component('comp1').physics('emw').prop('components').set('components', 'outofplane');
model.component('comp1').physics('emw').feature('wee1').set('epsilonr', {'epsrinf'; '0'; '0'; '0'; 'epsrinf'; '0'; '0'; '0'; 'epsrinf'});
model.component('comp1').physics('emw').feature('wee1').set('mur', {'murinf'; '0'; '0'; '0'; 'murinf'; '0'; '0'; '0'; 'murinf'});
model.component('comp1').physics('emw').feature('wee1').label('Wave Equation, Electric sca');
model.component('comp1').physics('emw').feature('wee2').set('epsilonr', {'epsilonb'; '0'; '0'; '0'; 'epsilonb'; '0'; '0'; '0'; 'epsilonb'});
model.component('comp1').physics('emw').feature('wee2').label('Wave Equation, Electric background');
model.component('comp1').physics('emw').feature('weak1').set('weakExpression', '-M1x*test(emw.curlEx)-M1y*test(emw.curlEy)');
model.component('comp1').physics('w').prop('Units').set('CustomSourceTermUnit', 'V^2/m^4');
model.component('comp1').physics('w').feature('wfeq1').set('weak', {'(test(M1x)*emw.curlEx+test(M1y)*emw.curlEy)*(omegam1*(i*alpha1*QNM_omega+omega01+omegam1))/nomega^2+(test(M1x)*emw.curlEy-test(M1y)*emw.curlEx)*(-i)*(QNM_omega*omegam1)/nomega^2-(test(M1x)*M1x+test(M1y)*M1y)*((1+alpha1^2)*QNM_omega^2-2*i*alpha1*QNM_omega*(omega01+omegam1)-(omega01+omegam1)^2)*murinf/nomega^2'; '0'});
model.component('comp1').physics('emw2').prop('components').set('components', 'outofplane');
model.component('comp1').physics('emw2').feature('wee1').set('epsilonr', {'epsrinf'; '0'; '0'; '0'; 'epsrinf'; '0'; '0'; '0'; 'epsrinf'});
model.component('comp1').physics('emw2').feature('wee1').set('mur', {'murinf'; '0'; '0'; '0'; 'murinf'; '0'; '0'; '0'; 'murinf'});
model.component('comp1').physics('emw2').feature('wee2').set('epsilonr', {'epsilonb'; '0'; '0'; '0'; 'epsilonb'; '0'; '0'; '0'; 'epsilonb'});
model.component('comp1').physics('emw2').feature('weak1').set('weakExpression', '-N1x*test(emw2.curlEx)-N1y*test(emw2.curlEy)');
model.component('comp1').physics('w2').prop('Units').set('CustomSourceTermUnit', 'V^2/m^4');
model.component('comp1').physics('w2').feature('wfeq1').set('weak', {'(test(N1x)*emw2.curlEx+test(N1y)*emw2.curlEy)*(omegam1*(i*alpha1*QNM_omega2+omega01+omegam1))/nomega^2+(test(N1x)*emw2.curlEy-test(N1y)*emw2.curlEx)*(i)*(QNM_omega2*omegam1)/nomega^2-(test(N1x)*N1x+test(N1y)*N1y)*((1+alpha1^2)*QNM_omega2^2-2*i*alpha1*QNM_omega2*(omega01+omegam1)-(omega01+omegam1)^2)*murinf/nomega^2'; '0'});

model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', '0.0097500');
model.component('comp1').mesh('mesh1').feature('size').set('hmin', '5.8500E-6');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmax', 'r/10');
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hmin', 5.85E-6);
model.component('comp1').mesh('mesh1').feature('ftri2').feature('size1').set('hminactive', false);
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hmax', 'r/2');
model.component('comp1').mesh('mesh1').feature('ftri3').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').run;

model.component('comp1').physics('emw').feature('wee1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee1').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('mur_mat', 'userdef');
model.component('comp1').physics('emw').feature('wee2').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw2').feature('wee1').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw2').feature('wee1').set('mur_mat', 'userdef');
model.component('comp1').physics('emw2').feature('wee1').set('sigma_mat', 'userdef');
model.component('comp1').physics('emw2').feature('wee2').set('epsilonr_mat', 'userdef');
model.component('comp1').physics('emw2').feature('wee2').set('mur_mat', 'userdef');
model.component('comp1').physics('emw2').feature('wee2').set('sigma_mat', 'userdef');

model.study.create('std1');
model.study('std1').create('eig', 'Eigenfrequency');
model.study('std1').feature('eig').set('activate', {'emw' 'on' 'w' 'on' 'emw2' 'off' 'w2' 'off' 'frame:spatial1' 'on'});
model.study.create('std2');
model.study('std2').create('eig', 'Eigenfrequency');
model.study('std2').feature('eig').set('activate', {'emw' 'off' 'w' 'off' 'emw2' 'on' 'w2' 'on' 'frame:spatial1' 'on'});

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('e1', 'Eigenvalue');
model.sol('sol1').feature('e1').create('d1', 'Direct');
model.sol.create('sol2');
model.sol('sol2').study('std2');
model.sol('sol2').attach('std2');
model.sol('sol2').create('st1', 'StudyStep');
model.sol('sol2').create('v1', 'Variables');
model.sol('sol2').create('e1', 'Eigenvalue');
model.sol('sol2').feature('e1').create('d1', 'Direct');

model.result.dataset.create('join1', 'Join');
model.result.dataset('join1').set('data', 'dset1');
model.result.dataset('join1').set('data2', 'dset2');
model.result.numerical.create('int1', 'IntSurface');
model.result.numerical('int1').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg2').create('surf1', 'Surface');

model.study('std1').feature('eig').set('neigs', 3);
model.study('std1').feature('eig').set('neigsactive', true);
model.study('std1').feature('eig').set('shift', 'freqg');
model.study('std2').feature('eig').set('neigs', 3);
model.study('std2').feature('eig').set('neigsactive', true);
model.study('std2').feature('eig').set('shift', 'freqg');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol1').feature('e1').set('neigs', 3);
model.sol('sol1').feature('e1').set('shift', 'freqg');
model.sol('sol1').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol1').runAll;
model.sol('sol2').attach('std2');
model.sol('sol2').feature('e1').set('transform', 'eigenfrequency');
model.sol('sol2').feature('e1').set('neigs', 3);
model.sol('sol2').feature('e1').set('shift', 'freqg');
model.sol('sol2').feature('e1').feature('aDef').set('complexfun', true);
model.sol('sol2').runAll;

model.result.dataset('join1').set('solutions', 'one');
model.result.dataset('join1').set('solnum', 2);
model.result.dataset('join1').set('solutions2', 'one');
model.result.dataset('join1').set('solnum2', 2);
model.result.dataset('join1').set('method', 'explicit');
model.result.numerical('int1').set('data', 'join1');
model.result.numerical('int1').set('solrepresentation', 'solnum');
model.result.numerical('int1').set('expr', {'(data1(emw.Ex)*data2(emw2.Dx)+data1(emw.Ey)*data2(emw2.Dy)+data1(emw.Ez)*data2(emw2.Dz))*data1(pml1.detInvT)' '-(data1(HLx)*data2(dwudwH_x)+data1(HLy)*data2(dwudwH_y)+data1(HLz)*data2(dwudwH_z))*data1(pml1.detInvT)*mu0_const' '(data1(emw.Ex)*data2(emw2.Dx)+data1(emw.Ey)*data2(emw2.Dy)+data1(emw.Ez)*data2(emw2.Dz))*data1(pml1.detInvT)-(data1(HLx)*data2(dwudwH_x)+data1(HLy)*data2(dwudwH_y)+data1(HLz)*data2(dwudwH_z))*data1(pml1.detInvT)*mu0_const' '-(data1(HLx)*data2(duwHx)+data1(HLy)*data2(duwHy)+data1(HLz)*data2(duwHz))*mu0_const*data1(pml1.detInvT)'});
model.result.numerical('int1').set('unit', {'N' 'N' 'N' 'N'});
model.result.numerical('int1').set('descr', {'QN_E' 'QN_H' 'QN=QN_E+QN_H' 'QN_H_m2'});
model.result('pg1').label('Electric Field (emw)');
model.result('pg1').set('looplevel', [6]);
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('surf1').label('Surface');
model.result('pg1').feature('surf1').set('expr', 'imag(emw.Ez/sqrt(7.2462E-19-1.1766E-19i))');
model.result('pg1').feature('surf1').set('descr', 'imag(emw.Ez/sqrt(7.2462E-19-1.1766E-19i))');
model.result('pg1').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg1').feature('surf1').set('smooth', 'internal');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg2').label('Electric Field (emw2)');
model.result('pg2').set('data', 'join1');
model.result('pg2').set('solrepresentation', 'solnum');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature('surf1').label('Surface');
model.result('pg2').feature('surf1').set('expr', 'real(data1(emw.Ez)*data2(emw2.Ez)/(6.5289E-19-2.3793E-20i))');
model.result('pg2').feature('surf1').set('unit', 'kg^2*m^2/(s^6*A^2)');
model.result('pg2').feature('surf1').set('descr', 'real(data1(emw.Ez)*data2(emw2.Ez)/(6.5289E-19-2.3793E-20i))');
model.result('pg2').feature('surf1').set('colortable', 'RainbowLight');
model.result('pg2').feature('surf1').set('resolution', 'normal');

out = model;
