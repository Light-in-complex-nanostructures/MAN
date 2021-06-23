function[cms, LH]= QNMnorm(cms, LH, iter)
% Find information on the electric-dipole in COMSOL model sheet
% We check that user defined the point dipole by a real dipole moment (and not by amplitude & direction).
model = mphload(cms.save_model);
dip_mom_str=model.physics('emw').feature(cms.electric_point_dipole_name).getStringArray('Iop');
cms.dipole_mom(1)=0;
cms.dipole_mom(2)=0;
cms.dipole_mom(3)=str2num(dip_mom_str); % Jz

% get the scattered field on the dipole
LH.E_source_scatt=LH.E_source(1:iter-LH.itercrash,:)-cms.E_source_bg(1:iter-LH.itercrash,:);
% calculate the normalization coefficient (cf OE 2013)
LH.normal_coeff=sqrt(1i*(LH.omega-LH.pole_estimate(1:(end-1)))./(dot(repmat(cms.dipole_mom,iter-LH.itercrash,1),LH.E_source_scatt,2)).');

%% PLOT the normalized field (tested component) on evaluation point to show the convergence of normalization process
LH.tested_field_scatt=LH.tested_field_tot-cms.tested_field_bg(1:iter-LH.itercrash);
LH.tested_field_normalized=LH.normal_coeff.*LH.tested_field_scatt./sqrt(cms.sym_factor);
%LH.normal_coeff=LH.normal_coeff.*sign(real(LH.tested_field_normalized));
LH.tested_field_normalized=LH.tested_field_normalized.*sign(real(LH.tested_field_normalized));
figure;
subplot(2,1,1)
plot(real(LH.tested_field_normalized),'mo-'); hold on; plot(real(LH.tested_field_normalized(1:3)),'bo-');
xlabel('Iteration number');ylabel(['Re($\tilde{' cms.tested_field_comp '}$)'],'Interpreter','LaTex');
title(['Re($\tilde{' cms.tested_field_comp '}$) at evaluation point'],'Interpreter','LaTex');
ax=gca;set(ax,'XTick',1:iter);clear ax
legend('intermediate estimates','initial triplet','Location','best');
subplot(2,1,2)
plot(imag(LH.tested_field_normalized),'mo-'); hold on; plot(imag(LH.tested_field_normalized(1:3)),'bo-');
xlabel('Iteration number');ylabel(['Im($\tilde{' cms.tested_field_comp '}$)'],'Interpreter','LaTex');
title(['Im($\tilde{' cms.tested_field_comp '}$) at evaluation point'],'Interpreter','LaTex');
ax=gca;set(ax,'XTick',1:iter);clear ax
fprintf('\n Normalized value of the QNM at evaluation point\n');
fprintf('\n\t value : %1.15e + %1.15e I \n',real(LH.tested_field_normalized(end)), imag(LH.tested_field_normalized(end)));

end