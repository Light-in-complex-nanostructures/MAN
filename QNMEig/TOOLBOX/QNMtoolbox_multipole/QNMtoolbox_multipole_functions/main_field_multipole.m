%% Convention physical constants
CT.c = 299792458; % Speed of light
CT.epsilon0 = 8.854187817e-12; % Vacuum permittivity
CT.mu0 = 4e-7*pi; % Vacuum permeability
%% User defined parameters 
model_path{1}='L:\QNM\9 released'; % Path of the COMSOL file
model_name{1}='QNM_Dolmen-sym.mph'; % Name of the COMSOL file
PEC=[1 0 0;]; % Perfect Electric Conductor Planes
PMC=[0 0 1;]; % Perfect Magnetic Conductor Planes
scatter  ='sel3';  % Domains occupied by the Scatter(can be either a series of numbers or a string)
matvar{1}='var12'; % Material of the resonator (Correspond to the 'Varibles' in the COMSOL used to define omegap1, gamma1, omega01, and epsiloninf)
% Variables for multipole decomposition
nmax=2; % Truncated order n for the VSWFs.
Ndes_th1=10; % Variable to control the number of points generated on the circumscribing sphere for integration. 
Ndes_th2=10;
Ndes_ph1=10;
Ndes_ph2=10;
r= 120e-9; % Radius of the circumscribing sphere
% Not used
sub='';
NG         =false; 
%% 1. Caluclate the eigen solution related variables
for isca=1:length(model_path)
    Sca0=modelPost_obj_multi2(model_path{isca}, model_name{isca}, r, nmax, Ndes_th1, Ndes_th2, Ndes_ph1, Ndes_ph2, scatter, matvar, sub);
    Sca0=extract_sol(Sca0);
    Sca0=extract_point_power(Sca0);
    Sca0=extract_point_sym(Sca0, PEC(isca,:), PMC(isca,:));    
    Sca0 = anm_bnm_eigen_QN(Sca0, CT, PEC(isca,:), PMC(isca,:));
    if(isca==1)
        Sca=Sca0;
    else
        Sca=modelPost_obj_multi2.mergeSca(Sca,Sca0,isca);
    end
end
Sca = file_management(Sca, model_path, model_name, PEC, PMC);
Sca = multipole_LZX(Sca, CT);


