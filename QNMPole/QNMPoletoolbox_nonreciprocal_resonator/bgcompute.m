function [E_source_bg, tested_field_bg]=bgcompute(cms, LH, iter)

model = mphload(cms.save_model); % we reload the model to refresh it
% Run the COMSOL simulation at 'omega' without the resonator to obtain the background field
fprintf('\tCOMSOL background calculations on pole frequency estimate\n');
info1=mphmodel(model.component('comp1').physics('emw'));

info1=fieldnames(info1);
[temp, Nodes]=find(strncmp(info1,"wee",3));
Nodes=length(Nodes);
selbgname=cell(1,Nodes); selbgent=cell(1,Nodes);
for cc=1:Nodes
    infobg=mphgetselection(model.component('comp1').physics('emw').feature(['wee',num2str(cc)]));
    selbgname{cc}=infobg.named;
    selbgent{cc} =infobg.entities;
    if(strcmp(['wee',num2str(cc)],cms.background_material))
        model.component('comp1').physics('emw').feature(['wee',num2str(cc)]).selection.all;
    elseif(cc~=1)
        model.component('comp1').physics('emw').feature(['wee',num2str(cc)]).active(false);
    end
end

try
% We create a 2nd study to store data for the background field calculations
model.study.create(cms.bg_field_study_name);
model.study(cms.bg_field_study_name).create ('freq', 'Frequency');
model.study(cms.bg_field_study_name).feature('freq').activate('emw', true);
model.study(cms.bg_field_study_name).feature('freq').set('plist', 'frequency_L_QNM');
model.study(cms.bg_field_study_name).feature('freq').set('punit', 'Hz');
catch
end

for cciter=4:iter-LH.itercrash
    omega_L_QNM=real(LH.pole_estimate(cciter));
    omega=LH.pole_estimate(cciter);
    model.param.set('omega_L_QNM', [num2str(omega_L_QNM,'%12.8e'),'[Hz]']);
    model.param.set(cms.omega_var_name,[num2str(omega,'%12.15e'),'[Hz]']);
    model.geom(cms.geometry_name).run();
    model.mesh(cms.mesh_name).run();
    % Run the simulation
    model.study(cms.bg_field_study_name).run(); % at this stage,the Maxwell eqs have been solved
    
    % get the field on the evaluation point
    [temp]=mphinterp(model,{['emw.' cms.tested_field_comp]},'dataset','dset2','coord',cms.eval_point,'Complexout','on');
    tested_field_bg(cciter)=temp;
    clear temp
    
    % get the E field on the dipole
    num_dipole_entity= model.physics('emw').feature(cms.electric_point_dipole_name).selection().entities(0);
    ddd = mpheval(model,{'Ex','Ey','Ez'},'dataset','dset2','edim','point','selection',num_dipole_entity,'Complexout','on');
    E_source_bg(cciter,:)=[ddd.d1(:),ddd.d2(:),ddd.d3(:)];
end
% Reset the resonator materials
for cc=1:Nodes
    if(cc~=1)
        if(isempty(selbgname{cc}))
            model.component('comp1').physics('emw').feature(['wee',num2str(cc)]).selection.set(selbgent{cc})
        else
            model.component('comp1').physics('emw').feature(['wee',num2str(cc)]).selection.named(selbgname{cc})
        end
        model.component('comp1').physics('emw').feature(['wee',num2str(cc)]).active(true);
    end
end
mphsave(model,cms.save_model);
end