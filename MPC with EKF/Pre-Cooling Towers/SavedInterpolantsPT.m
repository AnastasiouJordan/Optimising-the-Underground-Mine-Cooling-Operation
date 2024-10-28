clc
clear

% PT_meas  = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','PT','Range','C3:Q62756');
% SD_meas  = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','SD','Range','C3:L62756');
% Env_meas = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','Env','Range','C3:I62756');

PT_meas  = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','PT','Range','C11670:Q62756');
SD_meas  = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','SD','Range','C11670:L62756');
Env_meas = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','Env','Range','C11670:I62756');

% Required

F_outSD     = SD_meas{:,8};  % L/s, SD outlet flowrate
T_outSD     = SD_meas{:,9};  % oC,  SD outlet temperature
F_AntiSurge = PT_meas{:,3};  % L/s, Anti-surge flowrate
T_AntiSurge = PT_meas{:,4};  % oC,  Anti-surge temperature
F_outPT     = PT_meas{:,10}; % L/s, PT outlet flowrate
H_Env       = Env_meas{:,5}; % %,   Humidity of the environment
T_Env       = Env_meas{:,4}; % oC,  Dry bulb temp of the environment
T_PT        = PT_meas{:,9};  % oC,  PT outlet temperature, also calculated as v.T_PT
L_PTA       = PT_meas{:,5};  % %,   PT A level, also calculated as v.L_PTA
L_PTB       = PT_meas{:,6};  % %,   PT B level, also calculated as v.L_PTB


% Redundant
U_outPT     = PT_meas{:,8};  % %,   Valve position on the PT outlet
U_MP        = PT_meas{:,1};  % %,   Valve position at the mixing point

T = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','Env','Range','B11670:Q62756');
t = T{:,1};
t = t - 700020;

u.F_outSD     = griddedInterpolant(t, F_outSD);
u.T_outSD     = griddedInterpolant(t, T_outSD);
u.F_AntiSurge = griddedInterpolant(t, F_AntiSurge);
u.T_AntiSurge = griddedInterpolant(t, T_AntiSurge);
u.F_outPT     = griddedInterpolant(t, F_outPT);
u.H_Env       = griddedInterpolant(t, H_Env);
u.T_Env       = griddedInterpolant(t, T_Env);
u.T_PT        = griddedInterpolant(t, T_PT);
u.L_PTA       = griddedInterpolant(t, L_PTA );
u.L_PTB       = griddedInterpolant(t, L_PTB);
u.U_outPT     = griddedInterpolant(t, U_outPT);
u.U_MP        = griddedInterpolant(t, U_MP);

n.exogenousfields = {'F_outSD', 'T_outSD', 'F_AntiSurge',...
                     'T_AntiSurge', 'F_outPT', 'H_Env'...
                     'T_Env', 'T_PT', 'L_PTA', 'L_PTB'...
                     'U_outPT', 'U_MP'};

save SavedInterpolantsPT.mat u n t

clear all






