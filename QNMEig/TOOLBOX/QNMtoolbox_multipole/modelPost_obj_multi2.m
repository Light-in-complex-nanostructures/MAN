classdef modelPost_obj_multi2    
    properties
        model;
        % Variables
        in=struct('model_path',{},'model_name',{},'r',{},'nmax',{},'Ndes_th1',{},'Ndes_th2',{}, ...
            'Ndes_ph1',{},'Ndes_ph2',{},'epsln_b',{},'namb',{},'Nmaterial',{},'matvar',{},'scatter',{},'Nsca',{},'scamat',{},'scaeps',{},...
            'sub',{},'Nsub',{},'submat',{},'subeps',{});
        % Materials
        material=struct('npole',{},'gamma',{},'omega_p',{},'omega0',{},'epsln_inf',{},'entities',{});
        % Information about the input file
        inpart=struct('model_path',{},'model_name',{},'sols',{},'sole',{},'solnum',{},'Nfiles',{});
        % Symmetry planes
        in_sym=struct('PEC',{},'PMC',{});  
        % Solution information
        Sol=struct('solnum',{},'tb_freq',{},'QN',{});
        % Output multipole moments and coefficients computed with 
        % Eqs.(2-3)
        VSHe=struct('px_qnm',{}, 'py_qnm', {}, 'pz_qnm', {}, ...
            'mx_qnm',{}, 'my_qnm', {}, 'mz_qnm', {}, ...
            'EQ_qnm',{}, 'anm_qnm',{}, 'bnm_qnm',{});
        % Multipole moments in Cartesian coordinate computed with Eqs. (4-6) and those in  
        % ref.[Phys. Rev. B, 96, 035443 (2017)]
        CTSe= struct('pxyz_qnm_Lzx',{},'mxyz_qnm_Lzx',{},'txyz_qnm_Lzx',{}, ...
            'EQ_qnm',{},'MQ_qnm',{},'EO_qnm',{});
        % Struct to store the field distribution inside the scatter
        QNM;    
        
        QNM2;
        
    end
    %------------------------------------------------------
    methods
        function Y=modelPost_obj_multi2(model_path, model_name, r, nmax, Ndes_th1, Ndes_th2, Ndes_ph1, Ndes_ph2, scatter, matvar, sub)
            %  Read the input parameters assigned by the user.
            Y.in(1).model_path=model_path;
            Y.in(1).model_name=model_name;
            Y.in(1).r=r;
            Y.in(1).nmax=nmax;
            Y.in(1).Ndes_th1=Ndes_th1;
            Y.in(1).Ndes_ph1=Ndes_ph1;
            Y.in(1).Ndes_th2=Ndes_th2;
            Y.in(1).Ndes_ph2=Ndes_ph2;
            Y.in(1).scatter=scatter;
            Y.in(1).matvar=matvar;
            Y.in(1).sub=sub;
            
            modelfile=fullfile(Y.in(1).model_path,Y.in(1).model_name);
            Y.model = mphopen(modelfile);
            
            if (ischar(Y.in.scatter)==true)
                var= Y.model.component('mod1').selection(Y.in.scatter);
                info = mphgetselection(var);
                Y.in.scatter=info.entities;  % Domain numbers
            end
            Y.in(1).Nsca=length(Y.in.scatter);
            
            if (ischar(Y.in.sub)==true)
                if(strcmp(Y.in.sub,'')==false)
                    var= Y.model.component('mod1').selection(Y.in.sub);
                    info = mphgetselection(var);
                    Y.in.sub=info.entities;  % Domain numbers
                end
            end
            Y.in(1).Nsub=length(Y.in.sub);
               
        end
        %------------------------------------------------------
        function [Y]=extract_sol(Y)
            %  Extract, from the COMSOL file, the following parameters:
            % 1. The material parameters: omegap1, gamma1, omega01,
            % epsiloninf, omegap2, gamma2, omega02.
            % 2. The eigenfrequencies 'tb_freq'. (in COMSOL convention exp(i¦Øt))
            % 3. The modal normalizers 'QN'.(in COMSOL convention exp(i¦Øt))
            % 4. The background permittivity epsln_b.
            % 5. The permittivity of the domains scaeps.
            % For domains which are dispersive less scaeps gives
            % permittivity. (in COMSOL convention exp(i¦Øt))
            % For domains which are dispersive scaeps gives epsiloninf.
                
            model = Y.model;
            
            matvar=Y.in.matvar;
            Nmaterial =length(matvar); % Number of materials
            Y.in.Nmaterial=Nmaterial;
            Y.in.epsln_b=mphevaluate(model,'epsilonb');
            Y.in.namb=sqrt(Y.in.epsln_b);
            
            for imat=1:Nmaterial
                %    material=struct('npole',{},'gamma',{},'omega_p',{},'omega0',{},'epsln_inf',{},'entities',{});
                var=model.component('mod1').variable(matvar{imat});
                info = mphgetselection(var);
                try
                    Y.material(imat).entities=info.entities; % material domains
                catch
                    Y.material(imat).entities=Y.in.scatter; % All the scatterer has the same material
                end
                expr = mphgetexpressions(var);  % Expression of the material
                Y.material(imat).npole=(length(expr(:,1))-1)/3;  % Calculate the number of poles of the Drude model
                for ipole=1:Y.material(imat).npole
                    index = find(strcmp(expr(:,1), ['omegap',num2str(ipole)])); %
                    Y.material(imat).omega_p(ipole)=mphevaluate(model,expr{index,2});
                    
                    index = find(strcmp(expr(:,1), ['gamma',num2str(ipole)])); %
                    Y.material(imat).gamma(ipole)=mphevaluate(model,expr{index,2});
                    
                    index = find(strcmp(expr(:,1), ['omega0',num2str(ipole)])); %
                    Y.material(imat).omega0(ipole)=mphevaluate(model,expr{index,2});
                end
                index = find(strcmp(expr(:,1), 'epsiloninf')); %
                Y.material(imat).epsln_inf=mphevaluate(model,expr{index,2});
            end
            
            [Y.Sol(1).tb_freq]=mphglobal(model,{'freq'},'solnum','all','dataset','dset1','Complexout','on');
            [Y.Sol(1).QN]     =mphglobal(model,{'QN'},'solnum','all','dataset','dset1','Complexout','on');
            Y.Sol(1).solnum   =length(Y.Sol(1).tb_freq);
                     
            scamat=zeros(1,Y.in(1).Nsca);
            scaeps=scamat;
            submat=zeros(1,Y.in(1).Nsub);
            subeps=submat;
            % The material of different parts
            for isca=1:Y.in(1).Nsca
                [imat,epsilon_t]=idmaterial(Y, Y.in.scatter(isca));
                scamat(isca)=imat;
                if(imat==0)
                    scaeps(isca)=epsilon_t;  
                else
                    scaeps(isca) =Y.material(imat).epsln_inf; 
                end
            end
            
            for isca=1:Y.in(1).Nsub
                [imat,epsilon_t]=idmaterial(Y, Y.in.sub(isca));
                submat(isca)=imat;
                if(imat==0)
                    subeps(isca)=epsilon_t;  % The conjugate is due to comsol
                else
                    subeps(isca) =Y.material(imat).epsln_inf; 
                end
            end
            
            Y.in(1).scamat=scamat;
            Y.in(1).scaeps=scaeps;
            Y.in(1).submat=submat;
            Y.in(1).subeps=subeps;
            
        end
        %------------------------------------------------------
        function [ Y ] = extract_point_power(Y)
            % Extract the modal fields inside the resonator from the COMSOL file.
                   
            model=Y.model;
            scatter=Y.in.scatter;
            sub=Y.in.sub;
            Nsca=Y.in.Nsca;
            Nsub=Y.in.Nsub;
            
            Ntot=Y.in.Nsca+Y.in.Nsub;
            if(Y.in.Nsub>=1)
                if(Y.in.submat(1)==0)
                  Ntot=Y.in.Nsca;
                end
            end
            
            for isca=1:Ntot
                if(isca<=Nsca)
                    id=scatter(isca);
                else
                    id=sub(isca-length(scatter));
                end
                temp=mpheval(model, 'meshvol*pml1.detInvT', 'solnum', 1,'pattern','gauss','selection', id); % do not forget the pml1.detInvT in the PML region
                QNM(isca).mesh_vol=temp.d1;   % mesh volume
                QNM(isca).cord=temp.p;        % mesh coordinates
                temp=mpheval(model, 'Ex', 'solnum', 'all','pattern','gauss','selection', id,'Complexout','on');
                QNM(isca).Ex=temp.d1;         % Electric field x component
                temp=mpheval(model, 'Ey', 'solnum', 'all','pattern','gauss','selection', id,'Complexout','on');
                QNM(isca).Ey=temp.d1;         % Electric field y component
                temp=mpheval(model, 'Ez', 'solnum', 'all','pattern','gauss','selection', id,'Complexout','on');
                QNM(isca).Ez=temp.d1;         % Electric field y component
                
                clear temp;
                
                sz=size(QNM(isca).Ex);
                QNM(isca).Exyz=zeros(3,sz(1),sz(2));
                QNM(isca).Exyz(1,:,:)=QNM(isca).Ex;
                QNM(isca).Exyz(2,:,:)=QNM(isca).Ey;
                QNM(isca).Exyz(3,:,:)=QNM(isca).Ez;
               
            end
            Y.QNM=QNM;
        end
        %------------------------------------------------------
        function [Y] = extract_point_sym(Y, PEC, PMC)
            % Obtain the uncalculated modal fields using symmetry.
            
            model=Y.model;
            QNM=Y.QNM;
            Nsca=Y.in.Nsca;
            Nsub=Y.in.Nsub;
            solnum=Y.Sol.solnum;
            PEMC=PEC+PMC;
                        
            Ntot=Y.in.Nsca+Y.in.Nsub;
            if(Y.in.Nsub>=1)
                if(Y.in.submat(1)==0) % no need to consider the slab when it is non-disperisve
                  Ntot=Y.in.Nsca;
                end
            end
            
            for isca=1:Ntot
                cord=QNM(isca).cord;
                mesh_vol=QNM(isca).mesh_vol;
                
                for isym=1:3
                    if(PEMC(isym)==1)
                        n=zeros(1,3);
                        n(isym)=1;                      
                        cord2=[cord(1,:).*(-1)^n(1);cord(2,:).*(-1)^n(2);cord(3,:).*(-1)^n(3)];
                        cord=[cord cord2];
                        mesh_vol=[mesh_vol,mesh_vol];
                    end
                end
                
                [signE,signH,cord2] = modelPost_obj_multi2.symconvert(cord',3, PEC, PMC); 
                Ex  =QNM(isca).Ex;
                Ey  =QNM(isca).Ey;
                Ez  =QNM(isca).Ez;
                for isym=1:sum(PEC+PMC)
                    Ex =[Ex Ex];
                    Ey =[Ey Ey];
                    Ez =[Ez Ez];
                end
                
                for isol=1:solnum
                 Ex(isol,:)=Ex(isol,:).*signE{1};
                 Ey(isol,:)=Ey(isol,:).*signE{2};
                 Ez(isol,:)=Ez(isol,:).*signE{3};
                end                
                
                Y.QNM(isca).Ex=Ex;
                Y.QNM(isca).Ey=Ey;
                Y.QNM(isca).Ez=Ez;
                Y.QNM(isca).mesh_vol=mesh_vol;
                Y.QNM(isca).cord=cord;
                
                sz=size(Ex);
                Y.QNM(isca).Exyz=zeros(3,sz(1),sz(2));
                Y.QNM(isca).Exyz(1,:,:)=Ex;
                Y.QNM(isca).Exyz(2,:,:)=Ey;
                Y.QNM(isca).Exyz(3,:,:)=Ez;
            end
                           
        end
        %------------------------------------------------------
        function [Y] = file_management(Y, model_path, model_name, PEC, PMC)
            % File management
            Y.inpart(1).Nfiles=length(model_path);
            
            for ifile=1:Y.inpart.Nfiles
                Y.inpart(1).model_path{ifile}=model_path{ifile};
                Y.inpart(1).model_name{ifile}=model_name{ifile};
                Y.in_sym(1).PEC{ifile}=PEC(ifile,:);
                Y.in_sym(1).PMC{ifile}=PMC(ifile,:);
            end
            
            if (ifile==1)
                 Y.inpart.sols=1;
                 Y.inpart.sole=Y.Sol.solnum;
                 Y.inpart.solnum=Y.Sol.solnum;
            end
        end
        %------------------------------------------------------
        function[Y]=reopen(Y)
            modelfile=fullfile(Y.in(1).model_path,Y.in(1).model_name);
            Y.model=mphopen(modelfile);
        end
        %------------------------------------------------------
        function[Y]=samemodes(Y)
            % Find and remove the repeated modes in different files
            
            Nfiles=Y.inpart.Nfiles;
            PEC=Y.in_sym.PEC;
            PMC=Y.in_sym.PMC;
            tb_freq=Y.Sol.tb_freq;
            QN=Y.Sol.QN;
            QNM=Y.QNM;
            QNM2=Y.QNM2;
            inpart=Y.inpart;
            
            type=zeros(1,Nfiles);
            % symmetry type of the file       
            for ifile=1:Nfiles    
                type(ifile)=PEC{ifile}(1)*2^2+PEC{ifile}(2)*2^1+PEC{ifile}(3)*2^0;
                type(ifile)=PMC{ifile}(1)*2^5+PMC{ifile}(2)*2^4+PMC{ifile}(3)*2^3+type(ifile);
            end
            Ntb_freq=[];
            NQN=[];
            
            % find the repeated modes
            for itype=unique(type)
                ord=find(type==itype);
                cont=0;
                tfreq=[];
                for iord=ord
                    cont=cont+1;
                    sols=inpart.sols(iord);
                    sole=inpart.sole(iord);
                    repeat(sols:sole)=false;
                    for isol=sols:sole
                        if(any(abs((tfreq-tb_freq(isol))/tb_freq(isol))<1e-3))
                            repeat(isol)=true; % The eigenfrequency shows two times for a same symmetry
                        end
                    end
                    tfreq=[tfreq;tb_freq(sols:sole)];
                end
            end
            
            
            % remove the repeated modes
            NSolnum=0;
            for ifile=1:Nfiles
                sols=inpart.sols(ifile);
                sole=inpart.sole(ifile);
                Nsols(ifile)=1+NSolnum;                
                for isol=sols:sole
                    if(repeat(isol)==false)
                        NSolnum=NSolnum+1;
                        Ntb_freq(NSolnum)=tb_freq(isol);
                        NQN(NSolnum)=QN(isol);
                        NQNM.Ex(NSolnum,:)=QNM.Ex(isol,:);
                        NQNM.Ey(NSolnum,:)=QNM.Ey(isol,:);
                        NQNM.Ez(NSolnum,:)=QNM.Ez(isol,:);
                        NQNM.Exyz(:,NSolnum,:)=QNM.Exyz(:,isol,:);
                        NQNM2.Ex(NSolnum,:)=QNM2.Ex(isol,:);
                        NQNM2.Ey(NSolnum,:)=QNM2.Ey(isol,:);
                        NQNM2.Ez(NSolnum,:)=QNM2.Ez(isol,:);
                        NQNM2.Hx(NSolnum,:)=QNM2.Hx(isol,:);
                        NQNM2.Hy(NSolnum,:)=QNM2.Hy(isol,:);
                        NQNM2.Hz(NSolnum,:)=QNM2.Hz(isol,:);
                        NQNM2.Exyz(:,NSolnum,:)=QNM2.Exyz(:,isol,:);
                    end
                end
                Nsole(ifile)=NSolnum;         
            end
            
            % rewrite the variables in the module
            for ifile=1:Nfiles
                Y.inpart.sols(ifile)=Nsols(ifile);
                Y.inpart.sole(ifile)=Nsole(ifile);
                Y.inpart.solnum(ifile)=Nsole(ifile)-Nsols(ifile)+1;           
            end
            Y.Sol.tb_freq=transpose(Ntb_freq);
            Y.Sol.QN=NQN;
            Y.Sol.solnum=NSolnum;
            Y.QNM.Ex=NQNM.Ex;
            Y.QNM.Ey=NQNM.Ey;
            Y.QNM.Ez=NQNM.Ez;
            Y.QNM.Exyz=NQNM.Exyz;
            Y.QNM2.Ex=NQNM2.Ex;
            Y.QNM2.Ey=NQNM2.Ey;
            Y.QNM2.Ez=NQNM2.Ez;
            Y.QNM2.Hx=NQNM2.Hx;
            Y.QNM2.Hy=NQNM2.Hy;
            Y.QNM2.Hz=NQNM2.Hz;
            Y.QNM2.Exyz=NQNM2.Exyz;
            
                 
        end
        %------------------------------------------------------
        function [Y1]= anm_bnm_eigen_QN(Y1, CT, PEC, PMC)
            % Use Eq. (2) and Eq. (3) to calculate the a_nm and b_nm, and then further calculate the Cartesian multipole moments, p,  m, and  Qe.
            % The outputs are saved in VSHe. (in convention exp(-i¦Øt) which is different with COMSOL)
            
            % input
            Ndes_th1=Y1.in(1).Ndes_th1;
            Ndes_th2=Y1.in(1).Ndes_th2;
            Ndes_ph1=Y1.in(1).Ndes_ph1;
            Ndes_ph2=Y1.in(1).Ndes_ph2;
            r=Y1.in(1).r;
            solnum=Y1.Sol(1).solnum;
            tb_freq=Y1.Sol(1).tb_freq;
            model=Y1.model;
            namb=Y1.in(1).namb;
            nmax=Y1.in(1).nmax;
            QN=Y1.Sol(1).QN;
            
            % input for symmetry
            PEMC=PEC+PMC;
            
            % initialize
            anm_qnm  = zeros(nmax,2*nmax+1,solnum); bnm_qnm = anm_qnm;
            % multipole momentum in Cartesian coordinate
            px_qnm = zeros(1,solnum); py_qnm= px_qnm; pz_qnm= px_qnm;
            mx_qnm = zeros(1,solnum); my_qnm= mx_qnm; mz_qnm= mx_qnm;
            EQ_qnm = zeros(3,3,solnum);
                       
            % part 1: generate the points on the imaginary sphere
            [theta,dtheta]=retgauss(0,pi,Ndes_th1,Ndes_th2);
            [phi,dphi]=retgauss(0,2*pi,Ndes_ph1,Ndes_ph2);
            [T,P,R]=meshgrid(theta,phi,r); [dT,dP]=meshgrid(dtheta,dphi);
            X =R.*sin(T).*cos(P);
            Y =R.*sin(T).*sin(P);
            Z =R.*cos(T);
            
            cord=[X(:),Y(:),Z(:)];
            [signE,signH,cord2] = modelPost_obj_multi2.symconvert(cord,3,PEC,PMC);
            
            for isol=1:solnum
                % start to calculate the expansion coefficients for every solution
                % evaluate eigen fields (need to check)
                %  coord=[X(:),Y(:),Z(:)];
               
                [Esx,Esy,Esz]=mphinterp(model,{'emw.Ex','emw.Ey','emw.Ez'}, ...
                    'coord',cord2.','Complexout','on','dataset','dset1','solnum',isol);
                Esx = Esx/sqrt(QN(isol)).*signE{1};
                Esy = Esy/sqrt(QN(isol)).*signE{2};
                Esz = Esz/sqrt(QN(isol)).*signE{3};
                % QNM Projection on Mie vectors...
                n_omegam =conj(tb_freq(:,1)*2*pi); % the conj is due to comsol
                for n = 1:nmax
                    for m = -n:n
                        k0_qnm = namb*n_omegam(isol)/CT.c;
                        E0=1;
                        Emn = (abs(E0))/(2*sqrt(pi))*1i^(n-1+2*m)*sqrt((2*n+1)*factorial(n-m)/factorial(n+m));
                        
                        T=T(:); P=P(:); R=R(:); dT=dT(:); dP=dP(:); Esx=Esx(:); Esy=Esy(:); Esz=Esz(:);
                        
                        Er= Esx.*sin(T).*cos(P) + Esy.*sin(T).*sin(P)+Esz.*cos(T);
                        Et= Esx.*cos(T).*cos(P) + Esy.*cos(T).*sin(P)-Esz.*sin(T);
                        Ep=-Esx.*sin(P) + Esy.*cos(P);
                        
                        Er=Er(:); Et=Et(:); Ep=Ep(:);
                        
                        %Spherical Harmonic functions calculation (n,m)
                        
                        [Nr,Nt,Np] = Nnm_fixed(n, m, T, P, R, k0_qnm);
                        [Mr,Mt,Mp] = Mnm_fixed(n, m, T, P, R, k0_qnm);
                        
                        %coefficients of vector spherical harmonics
                        anm_qnm(n,nmax+1+m,isol) = 1/(Emn*k0_qnm^2)*...
                            sum(sin(T).*dT.*dP.*conj(Er.*Nr + Et.*Nt + Ep.*Np))/...
                            sum(sin(T).*dT.*dP.*(Nr.*conj(Nr) + Nt.*conj(Nt) + Np.*conj(Np)));
                        
                        bnm_qnm(n,nmax+1+m,isol) = 1/(Emn*k0_qnm^2)*...
                            sum(sin(T).*dT.*dP.*conj(Er.*Mr + Et.*Mt + Ep.*Mp))/...
                            sum(sin(T).*dT.*dP.*(Mr.*conj(Mr) + Mt.*conj(Mt) + Mp.*conj(Mp)));
                    end
                end
                
                [px_qnm(isol),py_qnm(isol),pz_qnm(isol),mx_qnm(isol),my_qnm(isol),mz_qnm(isol)]=dipole_test2(anm_qnm(:,:,isol),bnm_qnm(:,:,isol),nmax,k0_qnm,namb);
                [ EQ_qnm(1,1,isol), EQ_qnm(1,2,isol), EQ_qnm(1,3,isol), EQ_qnm(2,1,isol), EQ_qnm(2,2,isol), EQ_qnm(2,3,isol), EQ_qnm(3,1,isol), EQ_qnm(3,2,isol), EQ_qnm(3,3,isol)] = quadrupole_test( anm_qnm(:,:,isol), nmax, k0_qnm, namb );
                
            end
            
            Y1.VSHe(1).px_qnm=px_qnm;
            Y1.VSHe(1).py_qnm=py_qnm;
            Y1.VSHe(1).pz_qnm=pz_qnm;
            Y1.VSHe(1).mx_qnm=mx_qnm;
            Y1.VSHe(1).my_qnm=my_qnm;
            Y1.VSHe(1).mz_qnm=mz_qnm;
            Y1.VSHe(1).anm_qnm=anm_qnm;
            Y1.VSHe(1).bnm_qnm=bnm_qnm;
            
            Y1.VSHe(1).EQ_qnm=EQ_qnm;
            
            % electric dipole and magnetic dipole moments in cartesian coordinates
            % Ref. Metamaterials 5.2-3 (2011): 64-73.
            
        end
        %------------------------------------------------------
        function [ Y ] = multipole_LZX(Y, CT)
            % This subroutine calculates the multipoles through the method provided by
            % [1] ref. Optics letters 44 57 (2019) and
            % [2] ref. Physical Review B 96, 035443 (2017).
             % The outputs are all saved in CTSe. (in convention exp(-i¦Øt))
            
            % initialize
            model=Y.model;
            solnum=Y.Sol(1).solnum;
            scatter=Y.in(1).scatter;
            tb_freq=Y.Sol(1).tb_freq;
            epsln_b=Y.in(1).epsln_b;
            QN=Y.Sol(1).QN;
            mat=Y.in.scamat;
            eps=Y.in.scaeps;
            
            % intermadiate variables
            O_tensor=zeros(3,3,3,solnum);
            V_vector=zeros(3,solnum);
            QNM=Y.QNM;
            n_E_int  =zeros(3,solnum);
            n_Er_int =zeros(3,3,solnum);
            n_Err_int =zeros(3,3,3,solnum);
            
            %inline functions
            function y=no1(ixyz)
                if(ixyz==1)
                    y=2;
                elseif(ixyz==2)
                    y=3;
                elseif(ixyz==3)
                    y=1;
                end
            end
            
            function y=no2(ixyz)
                if(ixyz==1)
                    y=3;
                elseif(ixyz==2)
                    y=1;
                elseif(ixyz==3)
                    y=2;
                end
            end
            
            for isca=1:Y.in.Nsca
                
                pxyz_qnm_Lzx = zeros(3,solnum);
                mxyz_qnm_Lzx = zeros(3,solnum);
                txyz_qnm_Lzx = zeros(3,solnum);
                EQ_qnm = zeros(3,3,solnum);
                MQ_qnm = zeros(3,3,solnum);
                EO_qnm = zeros(3,3,3,solnum);
                
                n_omegam =  conj(tb_freq(:,1)*2*pi);  % The conj is due to comsol
                imat=mat(isca);
                if(imat==0) % the Constant material
                    n_epslnm= conj(eps(isca));  % The conjugate is due to comsol
                else
                    n_epslnm =Drude_epsln_multi( n_omegam, Y.material, imat); % The conjugate is die to exp(iwt)
                end
                
                % start to calculate the integration
                for ixyz1=1:3
                    if(solnum~=1)
                        n_E_int(ixyz1,:)=sum( bsxfun(@times,squeeze(QNM(isca).Exyz(ixyz1,:,:)),QNM(isca).mesh_vol),2);
                    else
                        t=reshape(QNM(isca).Exyz(ixyz1,:,:),[1,length(QNM(isca).mesh_vol)]);
                        n_E_int(ixyz1,:)=sum(t.*QNM(isca).mesh_vol);
                    end
                    
                    if(solnum~=1)
                        for ixyz2=1:3
                            n_Er_int(ixyz1, ixyz2, :)=sum( bsxfun(@times,squeeze(QNM(isca).Exyz(ixyz1,:,:)),QNM(isca).cord(ixyz2,:).*QNM(isca).mesh_vol),2);
                            for ixyz3=1:3
                                n_Err_int(ixyz1, ixyz2, ixyz3, :)= sum( bsxfun(@times,squeeze(QNM(isca).Exyz(ixyz1,:,:)),QNM(isca).cord(ixyz2,:).*QNM(isca).cord(ixyz3,:).*QNM(isca).mesh_vol),2);
                            end
                        end
                    else
                        for ixyz2=1:3
                            t=reshape(QNM(isca).Exyz(ixyz1,:,:),[1,length(QNM(isca).mesh_vol)]);
                            n_Er_int(ixyz1, ixyz2, :)=sum( t.*QNM(isca).cord(ixyz2,:).*QNM(isca).mesh_vol);
                            for ixyz3=1:3
                                n_Err_int(ixyz1, ixyz2, ixyz3, :)= sum( t.*QNM(isca).cord(ixyz2,:).*QNM(isca).cord(ixyz3,:).*QNM(isca).mesh_vol);
                            end
                        end
                    end
                end
                
                % normalize
                for isol=1:solnum
                    n_E_int(:,isol)=n_E_int(:,isol)/sqrt(QN(isol));
                    n_Er_int(:,:,isol)=n_Er_int(:,:,isol)/sqrt(QN(isol));
                    n_Err_int(:,:,:,isol)=n_Err_int(:,:,:,isol)/sqrt(QN(isol));
                end
                
                n_E_int  =conj(n_E_int);      % The conj is due to comsol
                n_Er_int =conj(n_Er_int);    % The conj is due to comsol
                n_Err_int=conj(n_Err_int);  % The conj is due to comsol
                
                fac0=transpose((n_epslnm(:)-epsln_b)*CT.epsilon0);
                fac_m=transpose(-1i*n_omegam/2);
                fac_t=transpose(n_omegam(:)*1i/10);
                fac_EQ=3.0;
                fac_MQ=transpose(n_omegam/3/1i);
                fac_V=1/5;
                
                for isol=1: solnum
                    % electric dipole
                    for ixyz=1:3
                        pxyz_qnm_Lzx(ixyz, isol)=fac0(isol)*n_E_int(ixyz,isol);
                    end
                    % magnetic dipole
                    mxyz_qnm_Lzx(1,isol)=fac_m(isol)*fac0(isol)*(n_Er_int(3,2,isol)-n_Er_int(2,3,isol));
                    mxyz_qnm_Lzx(2,isol)=fac_m(isol)*fac0(isol)*(n_Er_int(1,3,isol)-n_Er_int(3,1,isol));
                    mxyz_qnm_Lzx(3,isol)=fac_m(isol)*fac0(isol)*(n_Er_int(2,1,isol)-n_Er_int(1,2,isol));
                    % toroidal dipole
                    for ixyz=1:3
                        txyz_qnm_Lzx(ixyz,isol)=2*(n_Err_int(ixyz,1,1,isol)+n_Err_int(ixyz,2,2,isol)+n_Err_int(ixyz,3,3,isol));
                        txyz_qnm_Lzx(ixyz,isol)= txyz_qnm_Lzx(ixyz,isol) - ...
                            (n_Err_int(1,ixyz,1,isol)+n_Err_int(2,ixyz,2,isol)+n_Err_int(3,ixyz,3,isol));
                    end
                    txyz_qnm_Lzx(:,isol)=fac_t(isol)*fac0(isol)*txyz_qnm_Lzx(:,isol);
                    % electric quadrupole
                    EQ_com= n_Er_int(1,1,isol)+n_Er_int(2,2,isol)+n_Er_int(3,3,isol);
                    for ixyz1=1:3
                        for ixyz2=1:3
                            EQ_qnm(ixyz1,ixyz2,isol)=n_Er_int(ixyz1,ixyz2,isol)+n_Er_int(ixyz2,ixyz1,isol);
                        end
                        EQ_qnm(ixyz1,ixyz1,isol)= EQ_qnm(ixyz1,ixyz1,isol)-2/3*EQ_com;
                    end
                    EQ_qnm(:,:,isol)=EQ_qnm(:,:,isol)*fac_EQ*fac0(isol);
                    % magnetic quadrupole
                    for ixyz1=1:3
                        for ixyz2=1:3
                            MQ_qnm(ixyz1,ixyz2,isol)=n_Err_int(no2(ixyz1),no1(ixyz1),ixyz2,isol)-n_Err_int(no1(ixyz1),no2(ixyz1),ixyz2,isol);
                            MQ_qnm(ixyz1,ixyz2,isol)=MQ_qnm(ixyz1,ixyz2,isol)+ ...
                                n_Err_int(no2(ixyz2),no1(ixyz2),ixyz1,isol)-n_Err_int(no1(ixyz2),no2(ixyz2),ixyz1,isol);
                        end
                    end
                    MQ_qnm(:,:,isol)=fac_MQ(isol)*fac0(isol)*MQ_qnm(:,:,isol);
                    %electric octupole
                    for i0=1:3
                        for j0=1:3
                            for k0=1:3
                                O_tensor(i0,j0,k0,isol)=n_Err_int(i0,j0,k0,isol)+n_Err_int(j0,i0,k0,isol)+n_Err_int(k0,i0,j0,isol);
                            end
                        end
                    end
                    O_tensor(:,:,:,isol)=O_tensor(:,:,:,isol)*fac0(isol);
                    
                    for i0=1:3
                        V_vector(i0,isol)=2*(n_Err_int(1,1,i0,isol)+n_Err_int(2,2,i0,isol)+n_Err_int(3,3,i0,isol));
                        V_vector(i0,isol)=V_vector(i0,isol)+...
                            (n_Err_int(i0,1,1,isol)+n_Err_int(i0,2,2,isol)+n_Err_int(i0,3,3,isol));
                    end
                    V_vector(:,isol)=V_vector(:,isol)*fac0(isol)*fac_V;
                    
                    for i0=1:3
                        for j0=1:3
                            for k0=1:3
                                EO_qnm(i0,j0,k0,isol)=O_tensor(i0,j0,k0,isol);
                                EO_qnm(i0,j0,k0,isol)=EO_qnm(i0,j0,k0,isol)-kroneckerDelta_wt(i0,j0)*V_vector(k0);
                                EO_qnm(i0,j0,k0,isol)=EO_qnm(i0,j0,k0,isol)-kroneckerDelta_wt(i0,k0)*V_vector(j0);
                                EO_qnm(i0,j0,k0,isol)=EO_qnm(i0,j0,k0,isol)-kroneckerDelta_wt(j0,k0)*V_vector(i0);
                            end
                        end
                    end
                    
                end
                
                if(isca==1)
                    Y.CTSe(1).pxyz_qnm_Lzx = pxyz_qnm_Lzx;
                    Y.CTSe(1).mxyz_qnm_Lzx = mxyz_qnm_Lzx;
                    Y.CTSe(1).txyz_qnm_Lzx = txyz_qnm_Lzx;
                    Y.CTSe(1).EQ_qnm = EQ_qnm;
                    Y.CTSe(1).MQ_qnm = MQ_qnm;
                    Y.CTSe(1).EO_qnm = EO_qnm;
                else
                    Y.CTSe(1).pxyz_qnm_Lzx = pxyz_qnm_Lzx+Y.CTSe(1).pxyz_qnm_Lzx;
                    Y.CTSe(1).mxyz_qnm_Lzx = mxyz_qnm_Lzx+Y.CTSe(1).mxyz_qnm_Lzx;
                    Y.CTSe(1).txyz_qnm_Lzx = txyz_qnm_Lzx+Y.CTSe(1).txyz_qnm_Lzx;
                    Y.CTSe(1).EQ_qnm = EQ_qnm+Y.CTSe(1).EQ_qnm;
                    Y.CTSe(1).MQ_qnm = MQ_qnm+Y.CTSe(1).MQ_qnm;
                    Y.CTSe(1).EO_qnm = EO_qnm+Y.CTSe(1).EO_qnm;
                end
            end
            
        end
    end
    %------------------------------------------------------
    methods(Static)
        function [Y]=mergeSca(Y1,Y2,id)
            Y=Y1;
            Y.model=Y2.model;
            
            Y.Sol.solnum=Y1.Sol.solnum+Y2.Sol.solnum;
            Y.Sol.tb_freq=[Y1.Sol.tb_freq;Y2.Sol.tb_freq];
            Y.Sol.QN=[Y1.Sol.QN;Y2.Sol.QN];
            
            if(id==2)
                Y.inpart(1).solnum=Y1.Sol(1).solnum;
                Y.inpart(1).sols(1)=1;
                Y.inpart(1).sole(1)=Y1(1).Sol.solnum;
            end
            Y.inpart(1).solnum=[Y.inpart(1).solnum Y2.Sol.solnum];
            Y.inpart(1).sols(id)=Y.inpart(1).sole(id-1)+1;
            Y.inpart(1).sole(id)=Y.inpart(1).sols(id)+Y2.Sol.solnum-1;
            
            Ntot=Y.in.Nsca+Y.in.Nsub;
            if(Y.in.Nsub>=1)
                if(Y.in.submat(1)==0)
                  Ntot=Y.in.Nsca;
                end
            end
            
            for isca=1:Ntot
                Y.QNM(isca).Ex=[Y1.QNM(isca).Ex;Y2.QNM(isca).Ex];
                Y.QNM(isca).Ey=[Y1.QNM(isca).Ey;Y2.QNM(isca).Ey];
                Y.QNM(isca).Ez=[Y1.QNM(isca).Ez;Y2.QNM(isca).Ez];
                
                sz=size(Y.QNM(isca).Ex);
                Y.QNM(isca).Exyz=zeros(3,sz(1),sz(2));
                Y.QNM(isca).Exyz(1,:,:)=Y.QNM(isca).Ex;
                Y.QNM(isca).Exyz(2,:,:)=Y.QNM(isca).Ey;
                Y.QNM(isca).Exyz(3,:,:)=Y.QNM(isca).Ez;
            end
            
            for isca=1:Ntot
                dmesh=Y1.QNM(isca).mesh_vol-Y2.QNM(isca).mesh_vol;
                error=find(dmesh~=0);
                if(~isempty(error))
                    disp('error in merging files beacuse of different mesh')
                    return 
                end
            end
            
        end
        %------------------------------------------------------
        function [Y1] = extract_point_power_slice_sym( Y1, x, nx, y, ny, z, nz, dir, Nmx, Nmy, Nmz, dis, point)
             
            for ipart=1:Y1.inpart(1).Nfiles
                PEC=Y1.in_sym.PEC{ipart};
                PMC=Y1.in_sym.PMC{ipart};
                Y1.in.model_path=Y1.inpart(1).model_path{ipart};
                Y1.in.model_name=Y1.inpart(1).model_name{ipart};
                Y1=reopen(Y1);
                model=Y1.model;    
                domain=[];
                dmn=[];
                for imat=1:Y1.in(1).Nmaterial
                    domain=[domain Y1.material(imat).entities];
                    dmn=[dmn imat*ones(size(Y1.material(imat).entities))];
                end
                
                if(dir==3)
                 lx = linspace(nx, x, Nmx);
                 ly = linspace(ny, y, Nmy);
                 [X,Y] = meshgrid(lx,ly);
                 Z=ones(size(X))*dis;
                elseif(dir==1)
                 ly = linspace(ny, y, Nmy);
                 lz = linspace(nz, z, Nmz);
                 [Y,Z] = meshgrid(ly,lz);
                 X=ones(size(Y))*dis;
                elseif(dir==2)
                 lx = linspace(nx, x, Nmx);   
                 lz = linspace(nz, z, Nmz);    
                 [X,Z] = meshgrid(lx,lz);
                 Y=ones(size(Z))*dis;
                end
                 
                [sz1, sz2]=size(X);
                cx=reshape(X,[1,sz1*sz2]);
                cy=reshape(Y,[1,sz1*sz2]);
                cz=reshape(Z,[1,sz1*sz2]);
                coord=[cx',cy',cz'];
                          
                if(point==false)
                 [signE,signH,coord2] = modelPost_obj_multi2.symconvert(coord,3,PEC,PMC);
                else
                    coord2=coord;
                end
                            
                [Esx,Esy,Esz]=mphinterp(model,{'emw.Ex','emw.Ey','emw.Ez'}, ...
                    'coord',coord2.','Complexout','on','dataset','dset1','solnum','all');
                [Hsx,Hsy,Hsz]=mphinterp(model,{'emw.Hx','emw.Hy','emw.Hz'}, ...
                    'coord',coord2.','Complexout','on','dataset','dset1','solnum','all');
                
                if(point==false)
                    Esx=Esx.*signE{1};  Hsx=Hsx.*signH{1};
                    Esy=Esy.*signE{2};  Hsy=Hsy.*signH{2};
                    Esz=Esz.*signE{3};  Hsz=Hsz.*signH{3};
                end
                
                for isel=1:length(domain)
                    [DP1xA]=mphinterp(model,{'DP1x'}, ...
                        'coord',coord2.','Complexout','on','dataset','dset1','solnum',1,'selection',domain(isel));
                    if(isel==1)
                        id=zeros(size(DP1xA));
                    end
                    id(isnan(DP1xA)==false)=dmn(isel);
                end
                        
                if(ipart>1)
                    Y1.QNM2.Ex=[Y1.QNM2.Ex; Esx];
                    Y1.QNM2.Ey=[Y1.QNM2.Ey; Esy];
                    Y1.QNM2.Ez=[Y1.QNM2.Ez; Esz];
                    
                    Y1.QNM2.Hx=[Y1.QNM2.Hx; Hsx];
                    Y1.QNM2.Hy=[Y1.QNM2.Hy; Hsy];
                    Y1.QNM2.Hz=[Y1.QNM2.Hz; Hsz];
                else
                    Y1.QNM2.Ex=Esx;
                    Y1.QNM2.Ey=Esy;
                    Y1.QNM2.Ez=Esz;
                    
                    Y1.QNM2.Hx=Hsx;
                    Y1.QNM2.Hy=Hsy;
                    Y1.QNM2.Hz=Hsz;
                end
                Y1.QNM2.id=id;
                Y1.QNM2.cord=transpose(coord);
            end
            
            sz=size(Y1.QNM2.Ex);
            Y1.QNM2.Exyz=zeros(3,sz(1),sz(2));
            Y1.QNM2.Exyz(1,:,:)=Y1.QNM2.Ex;
            Y1.QNM2.Exyz(2,:,:)=Y1.QNM2.Ey;
            Y1.QNM2.Exyz(3,:,:)=Y1.QNM2.Ez;
            Y1.QNM2.Nmx=Nmx;
            Y1.QNM2.Nmy=Nmy;
            Y1.QNM2.Nmz=Nmz;
            Y1.QNM2.x=x; Y1.QNM2.nx=nx;
            Y1.QNM2.y=y; Y1.QNM2.ny=ny;
            Y1.QNM2.z=z; Y1.QNM2.nz=nz;
        end
        %------------------------------------------------------------------------------
        function [signE,signH,coord2] = symconvert(coord,dimension,PEC,PMC)
            
            PEMC=PEC+PMC;
            coord2=coord;
            for ixyz=1:dimension
                if(PEMC(ixyz)==1)
                    coord2(coord(:,ixyz)<0,ixyz)=-coord2(coord(:,ixyz)<0,ixyz);
                end
                signE{ixyz}=ones(1,length(coord2(:,1)));
                signH{ixyz}=ones(1,length(coord2(:,1)));
            end
            
            for ixyz=1:dimension
                if(PMC(ixyz)==1)
                    signE{ixyz}(coord(:,ixyz)<0)=signE{ixyz}(coord(:,ixyz)<0)*(-1);
                    for ixyz2=1:dimension
                        if(ixyz2~=ixyz)
                            signH{ixyz2}(coord(:,ixyz)<0)=signH{ixyz2}(coord(:,ixyz)<0)*(-1);
                        end
                    end
                elseif(PEC(ixyz)==1)
                    signH{ixyz}(coord(:,ixyz)<0)=signH{ixyz}(coord(:,ixyz)<0)*(-1);
                    for ixyz2=1:dimension
                        if(ixyz2~=ixyz)
                            signE{ixyz2}(coord(:,ixyz)<0)=signE{ixyz2}(coord(:,ixyz)<0)*(-1);
                        end
                    end
                end
            end      
        end
       %-----------------------------------------------------------
       
       
    end
end