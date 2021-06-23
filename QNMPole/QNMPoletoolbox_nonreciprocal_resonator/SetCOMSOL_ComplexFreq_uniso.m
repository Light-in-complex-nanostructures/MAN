function [model] = SetCOMSOL_ComplexFreq_uniso(model, omega_var_name, omega_L_QNM)
% This function automatically artificially allows COMSOL to perform calculations at complex
% frequencies by a trick which consists in modifying the Permittivity & Permeability distributions
% Details found in the Opt. Exp. paper, appendix 1
% The input variable "omega_L_QNM" should be real

study_name='std1';
model.param.set('omega_L_QNM', [num2str(omega_L_QNM,'%12.8e'),'[Hz]']);
model.param.set('frequency_L_QNM', 'omega_L_QNM/(2*pi)');

% get the number of EM nodes
info1=mphmodel(model.component('comp1').physics('emw'));
info1=fieldnames(info1);
[temp, Nodes]=find(strncmp(info1,"wee",3));
Nodes=length(Nodes);

% change the materials for each EM nods
for cc=1:Nodes
    info2 = mphgetproperties(model.component('comp1').physics('emw').feature(['wee',num2str(cc)]));
    % mu
    infomu=split(info2.mur,",");
    ePermeability=cell(1,length(infomu));
    for ccmu=1:length(infomu)
        if(not(strcmp( strtrim(infomu{ccmu}),'0')))
            infomu{ccmu}=strtrim(infomu{ccmu});
            ePermeability{ccmu}=['(',infomu{ccmu},')*',omega_var_name,'/omega_L_QNM'];
        else
            ePermeability{ccmu}='0';
        end
    end
    model.component('comp1').physics('emw').feature(['wee',num2str(cc)]').set('mur', ePermeability)
    
    %epsilon
    infoeps=split(info2.epsilonr,",");
    ePermittivity=cell(1,length(infoeps));
    for cceps=1:length(infoeps)
        if(not(strcmp( strtrim(infoeps{cceps}),'0')))
            infoeps{cceps}=strtrim(infoeps{cceps});
            ePermittivity{cceps}=['(',infoeps{cceps},')*',omega_var_name,'/omega_L_QNM'];
        else
            ePermittivity{cceps}='0';
        end
    end
    model.component('comp1').physics('emw').feature(['wee',num2str(cc)]').set('epsilonr', ePermittivity)
end

model.study(study_name).feature('freq').set('plist', 'frequency_L_QNM');

