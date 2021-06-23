function[signQNM]= QNMsign(cms, LH, RH, iter)

model = mphload(cms.save_model);
dip_mom_str=model.physics('emw').feature(cms.electric_point_dipole_name).getStringArray('Iop');
cms.dipole_mom(1)=0;
cms.dipole_mom(2)=0;
cms.dipole_mom(3)=str2num(dip_mom_str); % Jz

% get the scattered field on the dipole
E_source_scatt_temp=LH.E_source(1:iter-RH.itercrash,:)-cms.E_source_bg(1:iter-RH.itercrash,:);
% calculate the normalization coefficient (cf OE 2013)
normal_coeff=sqrt(1i*(RH.omega-RH.pole_estimate(1:(end-1)))./(dot(repmat(cms.dipole_mom,iter-RH.itercrash,1),E_source_scatt_temp,2)).');
% test QNM field when iter-RH.itercrash times of iteration are performed
lessit=LH.tested_field_scatt(iter-RH.itercrash)*normal_coeff(iter-RH.itercrash);
% test QNM field when iter-LH.itercrash times of iteration are performed
moreit=LH.tested_field_scatt(iter-LH.itercrash)*LH.normal_coeff(iter-LH.itercrash);
signQNM=sign(real(lessit*conj(moreit)));

end