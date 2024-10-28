clc
clear

% CD_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','CD','Range','C3:L62756');
% RP_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','RP','Range','C5:CG62758');
% Env_meas  = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','Env','Range','C3:I62756');

CD_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','CD','Range','C11670:L62756');
RP_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','RP','Range','C11672:CG62758');
Env_meas  = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','Env','Range','C11670:I62756');
PT_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','PT','Range','C11670:L62756');

% Required
L_CD    =  CD_meas{:,4}; % %,   Level in the CD
F_inRP  =  RP_meas{:,5}; % L/s, Inlet volmuetric flowrate to/through RPs
F_UG    =  CD_meas{:,9}; % L/s, Volumetric flowrate underground
F_outCD =  CD_meas{:,8}; % L/s, Outlet volumetric flowrate from the CD
T_CD    =  CD_meas{:,5}; % oC,  Outlet temperature from the CD
T_outRP =  RP_meas{:,7}; % oC,  Outlet temperature from the RPs
T_Env   = Env_meas{:,4}; % oC,  Dry bulb temp of the environment
F_outPT = PT_meas{:,10}; % L/s, Outlet volumetric flowrate from the PTs
T_outPT =  PT_meas{:,8}; % oC,  Outlet temperature from the PTs
L_PTA   =  PT_meas{:,5}; % %,   Level in the PT A basin
L_PTB   =  PT_meas{:,6}; % %,   Level in the PT B basin

% Redundant
T_inCD  = CD_meas{:,3};  % oC,  Inlet temperature to the CD
U_outCD = CD_meas{:,7};  % %,   Valve position on CD outlet
T_UG    = CD_meas{:,10}; % oC,  Temperature of flow underground

T = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','CD','Range','B11670:Q62756');
t = T{:,1};
t = t - 700020;

u.L_CD    = griddedInterpolant(t, L_CD);
u.F_inRP  = griddedInterpolant(t, F_inRP);
u.F_UG    = griddedInterpolant(t, F_UG);
u.F_outCD = griddedInterpolant(t, F_outCD);
u.T_CD    = griddedInterpolant(t, T_CD);
u.T_outRP = griddedInterpolant(t, T_outRP);
u.T_Env   = griddedInterpolant(t, T_Env);
u.F_outPT = griddedInterpolant(t, F_outPT);
u.T_outPT = griddedInterpolant(t, T_outPT);

u.T_inCD  = griddedInterpolant(t, T_inCD);
u.U_outCD = griddedInterpolant(t, U_outCD);
u.T_UG    = griddedInterpolant(t, T_UG);

F_Ice     = max((u.F_outPT(t)-u.F_UG(t)),0); % Flowrate to ice plants

u.F_Ice   = griddedInterpolant(t, F_Ice);

u.L_PTA   = griddedInterpolant(t, L_PTA);
u.L_PTB   = griddedInterpolant(t, L_PTB);

n.exogenousfields = {'L_CD', 'F_inRP', 'F_UG',...
                     'F_outCD', 'T_CD', 'T_outRP','T_Env',...
                     'F_outPT','T_outPT','T_inCD', 'U_outCD',... 
                     'T_UG','F_Ice','L_PTA', 'L_PTB'};

save SavedInterpolantsCD.mat u n t

clear all