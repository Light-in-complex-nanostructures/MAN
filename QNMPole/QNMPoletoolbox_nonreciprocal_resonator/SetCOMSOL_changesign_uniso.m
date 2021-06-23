function [model] = SetCOMSOL_changesign_uniso(model)
% This function automatically artificially allows COMSOL to perform calculations at complex
% frequencies by a trick which consists in modifying the Permittivity & Permeability distributions
% Details found in the Opt. Exp. paper, appendix 1
% The input variable "omega_L_QNM" should be real
% get the number of EM nodes
info1=mphmodel(model.component('comp1').physics('emw'));
info1=fieldnames(info1);
[temp, Nodes]=find(strncmp(info1,"wee",3));
Nodes=length(Nodes);
ng=[2 3 4 7 8];
% change the materials for each EM nods
for cc=1:Nodes 
    info2 = mphgetproperties(model.component('comp1').physics('emw').feature(['wee',num2str(cc)]));
    % mu   
    infomu=split(info2.mur,",");

    ePermeability=infomu;
    for ccmu=1:length(infomu)
        if(not(strcmp( strtrim(infomu{ccmu}),'0')))
            if(any(ng == ccmu))
                infomu{ccmu}=strtrim(infomu{ccmu});
                ePermeability{ccmu}=['(',infomu{ccmu},')*','(-1)'];
            end
        end
    end
    model.component('comp1').physics('emw').feature(['wee',num2str(cc)]').set('mur', ePermeability)
    
    %epsilon
    infoeps=split(info2.epsilonr,",");

    ePermittivity=infoeps;
    for cceps=1:length(infoeps)
        if(not(strcmp( strtrim(infoeps{cceps}),'0')))
            if(any(ng == cceps))
                infoeps{cceps}=strtrim(infoeps{cceps});
                ePermittivity{cceps}=['(',infoeps{cceps},')*','(-1)'];
            end
        end
    end
    model.component('comp1').physics('emw').feature(['wee',num2str(cc)]').set('epsilonr', ePermittivity)
end

