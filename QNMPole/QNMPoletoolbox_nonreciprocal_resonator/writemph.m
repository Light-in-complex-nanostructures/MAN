function []=writemph(cms, LH, RH, signQNM)

model = mphload(cms.save_model); % we reload the model to refresh it
model.result.param.set('omega_N',num2str(LH.pole_estimate(end),'%12.15e'));% send normalisation frequency back to save_model
model.result.param.set('N_coef',num2str(LH.normal_coeff(end),'%12.15e'));% send normalisation coefficient back to save_model
model.result.param.set('N_coefR',num2str(RH.normal_coeff(end),'%12.15e'));% send normalisation coefficient back to save_model
model.result.param.set('signQNM',num2str(signQNM,'%12.15e'));% send normalisation coefficient back to save_model
model.result.param.set('sym_factor',num2str(cms.sym_factor));% send symmetry factor back to save_model (not taken into account in N_coef !)


model.result.dataset.create('join1', 'Join');
model.result.dataset('join1').set('data', 'dset1');
model.result.dataset('join1').set('data2', 'dset3');
model.result.dataset('join1').set('method', 'explicit');
model.result.dataset('join1').set('solutions', 'one');
model.result.dataset('join1').set('solutions2', 'one');
model.result('pg1').feature('surf1').set('expr', 'data1(emw.Ez)*data2(emw.Ez)*N_coef*N_coefR*signQNM');
model.result('pg1').set('data', 'join1');

mphsave(model,cms.save_model)



end