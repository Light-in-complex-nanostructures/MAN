function [idmat,epsilon]=idmaterial(Sca, domain)
fd=[];
%----identify the material of the substrate-----
for imat=1:Sca.in.Nmaterial
    fd=find(Sca.material(imat).entities==domain);
    if(length(fd)==1)
        idmat=imat;
        break;
    end
end

if( isempty(fd))
    idmat=0;  % constant material  
    temp=mpheval(Sca.model, 'emw.epsilonrxx', 'solnum', 'all','pattern','gauss','selection', domain,'Complexout','on');
    epsilon=temp.d1(1);
else
    epsilon=Sca.material(idmat).epsln_inf;
end


