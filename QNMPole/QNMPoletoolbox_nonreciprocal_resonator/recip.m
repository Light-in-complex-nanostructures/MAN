function [RH]=recip(cms, LH, iter)

model = mphload(cms.save_model); % we reload the model to refresh it
% Run the COMSOL simulation at 'omega' without the resonator to obtain the background field
fprintf('\reciprocal mode estimate\n');
info1=mphmodel(model.component('comp1').physics('emw'));

% change the sign of the non-diagonal terms of the permittivity
[model] = SetCOMSOL_changesign_uniso(model);
try
    % We create a 3nd study to store data for the background field calculations
    model.study.create(cms.res_field_study_name);
    model.study(cms.res_field_study_name).create ('freq', 'Frequency');
    model.study(cms.res_field_study_name).feature('freq').activate('emw', true);
    model.study(cms.res_field_study_name).feature('freq').set('plist', 'frequency_L_QNM');
    model.study(cms.res_field_study_name).feature('freq').set('punit', 'Hz');
    model.sol('sol3').feature('s1').feature('dDef').active(true);
catch
end
RH.itercrash=LH.itercrash;
% recompute the fields for the reciprocal system
for cciter=1:iter-LH.itercrash
    omega_L_QNM=real(LH.pole_estimate(cciter));
    omega=LH.pole_estimate(cciter);
    model.param.set('omega_L_QNM', [num2str(omega_L_QNM,'%12.8e'),'[Hz]']);
    model.param.set(cms.omega_var_name,[num2str(omega,'%12.15e'),'[Hz]']);
    model.geom(cms.geometry_name).run();
    model.mesh(cms.mesh_name).run();
    % Run the simulation
    try
        model.study(cms.res_field_study_name).run(); % at this stage,the Maxwell eqs have been solved
    catch
        RH.itercrash=iter-cciter+1;
        %   RH.E_source(iter-LH.itercrash,:)=[0,0,0];
        %   RH.tested_field_tot(iter-LH.itercrash)=0;
        RH.omega=LH.pole_estimate(iter-RH.itercrash+1);
        break;
    end
    %We retrieve electric point dipole position
    num_dipole_entity= model.physics('emw').feature(cms.electric_point_dipole_name).selection().entities(0);%COMSOL point entity the dipole is attached to
    ddd = mpheval(model,{'Ex','Ey','Ez'},'edim','point','selection',num_dipole_entity,'Complexout','on','dataset','dset3');
    RH.E_source(cciter,:)=[ddd.d1(:),ddd.d2(:),ddd.d3(:)];
    
    [temp]=mphinterp(model,{['emw.' cms.tested_field_comp]},'coord',cms.eval_point,'Complexout','on','dataset','dset3');
    RH.tested_field_tot(cciter)=temp; clear t_field
    RH.pole_estimate(cciter)=LH.pole_estimate(cciter);
    RH.pole_estimate(cciter+1)=LH.pole_estimate(cciter+1);
    RH.omega=LH.pole_estimate(cciter+1);
    mphsave(model,cms.save_model)
    cciter
end

% change the sign back
model = mphload(cms.save_model); % we reload the model to refresh it
[model] = SetCOMSOL_changesign_uniso(model);
mphsave(model,cms.save_model)

end