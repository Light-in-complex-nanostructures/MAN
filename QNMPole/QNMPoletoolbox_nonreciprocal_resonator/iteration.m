function [E_source, tested_field_set, tested_field_tot, pole_estimate, time_cal, itercrash, omega_set, omega]=iteration( ...
                      model, omega, omega_var_name, geometry_name, mesh_name,tested_field_tot, h1,...
                      total_field_study_name, electric_point_dipole_name, tested_field_comp, eval_point, iter, ...
                      E_source, tested_field_set, pole_estimate, time_cal, omega_set, delta, save_model, QNM_ite)


%%
SetCOMSOL_ComplexFreq_uniso(model,omega_var_name,real(omega));
% Run the COMSOL simulation at 'omega'
fprintf('\tCOMSOL field calculations on pole frequency estimate\n');
% set the frequency of COMSOL computation
model.param.set(omega_var_name,[num2str(omega,'%12.15e'),'[Hz]']);
model.geom(geometry_name).run();
model.mesh(mesh_name).run();
itercrash=0;
% Run the simulation
if ~exist('debug_mode','var') % Matlab catches COMSOL exception to avoid program crash
    try
        model.sol('sol1').feature('s1').feature('dDef').active(true);
        model.study(total_field_study_name).run(); % at this stage,the Maxwell eqs have been solved
    catch
        fprintf('\n Iteration %s did not converge on users criteria: QNM pole reached (possibly)\n',num2str(iter));
        fprintf('\n\t Pole : %1.15e + %1.15e I \n',real(omega), imag(omega));
        itercrash=1;
        toc
        figure(h1);
        subplot(2,2,2);hold on;plot(iter,real(pole_estimate(end)),'ro','LineWidth',2);
        subplot(2,2,4);hold on;plot(iter,imag(pole_estimate(end)),'ro','LineWidth',2);
        return
    end
else
    model.study(total_field_study_name).run(); % at this stage,the Maxwell eqs have been solved
end
% get the field on the point at which the field is calculated for the iterative pole search
[temp]=mphinterp(model,{['emw.' tested_field_comp]},'coord',eval_point,'Complexout','on');
tested_field_tot(iter)=temp; clear t_field
clear temp

%We retrieve electric point dipole position
num_dipole_entity= model.physics('emw').feature(electric_point_dipole_name).selection().entities(0);%COMSOL point entity the dipole is attached to
% get the field on the dipole
ddd = mpheval(model,{'Ex','Ey','Ez'},'edim','point','selection',num_dipole_entity,'Complexout','on');
E_source(iter,:)=[ddd.d1(:),ddd.d2(:),ddd.d3(:)];
clear dip_mom_str
%% Compute new set of (3) interpolating frequencies
% generate the new frequency from triplet (see OE by Qiang Bai)
tested_field_set(end)=tested_field_tot(iter);% we add the newly calculated field in the set
[omega_set,tested_field_set]=omega_generation_anis(omega_set,tested_field_set,delta);
omega=omega_set(end);% the last value of the new omega set is the next to be calculated (new pole estimation)
pole_estimate(iter+1)=omega;  % save the frequency at each calculation
time_cal(iter)=toc;
%% Plot field value for each iteration (|Field| must diverge when we reach the pole)
figure(h1);
if iter>3
    % PLOT the field divergence
    subplot(1,2,1)
    semilogy(1:iter,1./abs(tested_field_tot),'-mo',1:3,1./abs(tested_field_tot(1:3)),'bo-');
    xlabel('Iteration number');ylabel(['1/|' tested_field_comp '| at evaluation point']);title('Field divergence as one approaches the pole');
    ax=gca;set(ax,'XTick',1:iter);clear ax
    legend('intermediate estimates','initial triplet','Location','southwest');
    % PLOT the pole frequency estimations convergence
    subplot(2,2,2)
    plot(real(pole_estimate),'mo-'); hold on; plot(real(pole_estimate(1:3)),'bo-'); hold off ;
    ax=gca;set(ax,'XTick',1:iter);clear ax
    title('Pole real part (rad/s)');
    subplot(2,2,4)
    plot(imag(pole_estimate),'mo-'); hold on; plot(imag(pole_estimate(1:3)),'bo-'); hold off ;
    ax=gca;set(ax,'XTick',1:iter);clear ax
    xlabel('Iteration number');title('Pole imaginary part (rad/s)');
    drawnow
else
    % PLOT the field divergence
    subplot(1,2,1)
    semilogy(1:iter,1./abs(tested_field_tot),'bo-');
    xlabel('Iteration number');ylabel(['1/|' tested_field_comp '| at evaluation point']);title('Field divergence as one approaches the pole');
    ax=gca;set(ax,'XTick',1:iter);clear ax
    %legend('initial triplet','Location','southwest');
    % PLOT the pole frequency estimations convergence
    subplot(2,2,2)
    plot(real(pole_estimate),'bo-');
    ax=gca;set(ax,'XTick',1:iter);clear ax
    title('Pole real part (rad/s)');
    %legend('intermediate estimates','initial triplet','pole','Location','best');
    subplot(2,2,4)
    plot(imag(pole_estimate),'bo-');
    ax=gca;set(ax,'XTick',1:iter);clear ax
    xlabel('Iteration number');title('Pole imaginary part (rad/s)');
    drawnow
end

%% we save the file for every iteration 
mphsave(model,save_model);% We save the model with all the data
if(iter==QNM_ite)
    subplot(2,2,2);hold on;plot(iter,real(pole_estimate(end)),'ro','LineWidth',2);
    subplot(2,2,4);hold on;plot(iter,imag(pole_estimate(end)),'ro','LineWidth',2);
end
end

