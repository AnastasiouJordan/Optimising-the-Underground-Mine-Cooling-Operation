clc
clear

T = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','CD','Range','B11670:Q62756');
t = T{:,1};
t = t - 700020;

% Read in data
SD_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','SD','Range','C11670:L62756');
PT_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','PT','Range','C11670:Q62756');
RP_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','RP','Range','C11672:CG62758');
CD_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','CD','Range','C11670:L62756');
Env_meas  = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','Env','Range','C11670:I62756');

% Square Dam
L_SD    =  SD_meas{:,6};  % %,   Level in the SD (also calculated)
F_inSD  =  SD_meas{:,5};  % L/s, Inlet volmuetric flowrate to the SD
F_outSD =  SD_meas{:,8};  % L/s, Outlet volumetric flowrate from the CD
T_outSD =  SD_meas{:,9};  % oC, Temperature of the water leaving the SD
U_outSD =  SD_meas{:,10}; % %,  Valve position on SD outlet

% Pre-cooling Towers
F_AntiSurge = PT_meas{:,3};  % L/s, Anti-surge flowrate
T_AntiSurge = PT_meas{:,4};  % oC,  Anti-surge temperature
F_outPT     = PT_meas{:,10}; % L/s, PT outlet flowrate
H_Env       = Env_meas{:,5}; % %,   Humidity of the environment
T_Env       = Env_meas{:,4}; % oC,  Dry bulb temp of the environment
T_PT        = PT_meas{:,9};  % oC,  PT outlet temperature, also calculated as v.T_PT
L_PTA       = PT_meas{:,5};  % %,   PT A level, also calculated as v.L_PTA
L_PTB       = PT_meas{:,6};  % %,   PT B level, also calculated as v.L_PTB
U_outPT     = PT_meas{:,8};  % %,   Valve position on the PT outlet
U_MP        = PT_meas{:,1};  % %,   Valve position at the mixing point

% Fridge Plants
F_inRP   =  RP_meas{:,5}; % L/s, Total inlet volumetric flowrate to RPs
T_outRP  =  RP_meas{:,7}; % oC,  Outlet temperature from the RPs

% Chill Dam
L_CD    =  CD_meas{:,4}; % %,   Level in the CD
T_CD    =  CD_meas{:,5}; % oC,  Outlet temperature from the CD
T_inCD  =  CD_meas{:,3}; % oC,  Inlet temperature to the CD
U_outCD =  CD_meas{:,7}; % %,   Valve position on CD outlet
F_outCD =  CD_meas{:,8}; % L/s, Outlet volumetric flowrate from the CD

% Underground
F_UG    =  CD_meas{:,9};  % L/s, Volumetric flowrate underground
T_UG    =  CD_meas{:,10}; % oC,  Temperature of flow underground

% Environment
T_Env   = Env_meas{:,4}; % oC,  Dry bulb temp of the environment

% Time series data
% Square Dam
u.L_SD    = griddedInterpolant(t, L_SD);
u.F_inSD  = griddedInterpolant(t, F_inSD);
u.F_outSD = griddedInterpolant(t, F_outSD);
u.T_outSD = griddedInterpolant(t, T_outSD);
u.U_outSD = griddedInterpolant(t, U_outSD);

% Pre-cooling Towers
u.F_AntiSurge = griddedInterpolant(t, F_AntiSurge);
u.T_AntiSurge = griddedInterpolant(t, T_AntiSurge);
u.F_outPT     = griddedInterpolant(t, F_outPT);
u.H_Env       = griddedInterpolant(t, H_Env);
u.T_PT        = griddedInterpolant(t, T_PT);
u.L_PTA       = griddedInterpolant(t, L_PTA );
u.L_PTB       = griddedInterpolant(t, L_PTB);
u.U_outPT     = griddedInterpolant(t, U_outPT);
u.U_MP        = griddedInterpolant(t, U_MP);


% Fridge Plants
u.F_inRPtot = griddedInterpolant(t, F_inRP);
u.T_amb     = griddedInterpolant(t, T_Env);
u.T_outPT   = griddedInterpolant(t, T_PT); 
u.T_outCD   = griddedInterpolant(t, T_CD);
% Per fridge plant, j
% On/off statuses
for j = 1:5
    k = (15*j) + 1;
    sj(:,j) = RP_meas{:, k};
end
sj(sj <= 0.1) = 0;
sj(sj > 0.1)  = 1;
u.s = griddedInterpolant(t, sj); 
% Inlet volumetric flowrates
for j = 1:5
    k = (15*j) - 2;
    F_inRPj(:,j) = RP_meas{:, k};
end
u.F_inRP = griddedInterpolant(t, F_inRPj); 
% Outlet volumetric flowrates
for j = 1:5
    k = (15*j);
    F_outRPj(:,j) = RP_meas{:, k}; 
end
u.F_outRP = griddedInterpolant(t, F_outRPj);
% Inlet temperatures (oC)
for j = 1:5
    k = (15*j) - 4;
    T_inRPj(:,j) = RP_meas{:, k};
end
u.T_inRP = griddedInterpolant(t, T_inRPj); 
% Outlet temperatures
for j = 1:5
    k = (15*j) - 1;
    T_RPj(:,j) = RP_meas{:, k};
end
u.T_RP = griddedInterpolant(t, T_RPj); 
% Refrigerant temperatures (liquid)
for j = 1:5
    k = (15*j) + 8;
    T_refrj(:,j) = RP_meas{:, k};
end
u.T_refr = griddedInterpolant(t, T_refrj); 


% Chill Dam
u.L_CD    = griddedInterpolant(t, L_CD);
u.F_outCD = griddedInterpolant(t, F_outCD);
u.T_CD    = griddedInterpolant(t, T_CD);
u.T_outRP = griddedInterpolant(t, T_outRP);
u.T_outPT = griddedInterpolant(t, T_PT);
u.T_inCD  = griddedInterpolant(t, T_inCD);
u.U_outCD = griddedInterpolant(t, U_outCD);

u.T_UG    = griddedInterpolant(t, T_UG);
u.F_UG    = griddedInterpolant(t, F_UG);
u.T_Env   = griddedInterpolant(t, T_Env);

F_Ice     = max((u.F_outPT(t)-u.F_UG(t)),0); % Flowrate to ice plants
u.F_Ice   = griddedInterpolant(t, F_Ice);

save SavedInterpolants.mat u t

clear all