function [model] = SetCOMSOL_ComplexFreq(model,omega_var_name,omega_L_QNM)
% This function automatically artificially allows COMSOL to perform calculations at complex  
% frequencies by a trick which consists in modifying the Permittivity & Permeability distributions 
% Details found in the Opt. Exp. paper, appendix 1
% The input variable "omega_L_QNM" should be real
study_name='std1';

model.param.set('omega_L_QNM', num2str(omega_L_QNM,'%12.8e'));
model.param.set('frequency_L_QNM', 'omega_L_QNM/(2*pi)');

Mats=model.material.tags();
for i=1:numel(Mats)
    additionnal=strcat(')*',omega_var_name,'/omega_L_QNM');
    Permittivity=char(model.material(Mats(i)).propertyGroup('def').getString('relpermittivity'));
    model.material(Mats(i)).propertyGroup('def').set('relpermittivity', {strcat('(',Permittivity,additionnal)});
    Permeability=char(model.material(Mats(i)).propertyGroup('def').getString('relpermeability'));
    model.material(Mats(i)).propertyGroup('def').set('relpermeability', {strcat('(',Permeability,additionnal)});
end
model.study(study_name).feature('freq').set('plist', 'frequency_L_QNM');
end

