%% Combined System: Square Dam, Pre-cooling Towers, Fridge Plants and Chill Dam
%  Jordan Anastasiou, 2023-12
%  This code is for the entire system, including
%  MPC and optimisation for the scheduling of
%  the five fridge plants according to Eskom's
%  time-of-day tariffs.

clc
clear
clf

%% Define exogeneous inputs
load SavedInterpolants.mat

% line = zeros(1, length(t));
% colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
% yyaxis left
% plot(t/86400, u.F_outCD(t),'color',[0 0.4470 0.7410], 'LineWidth',1); 
% hold on
% plot(t/86400, line', '-k', 'LineWidth', 1.5);
% ylabel('F_o_u_t_C_D (L/s)'); xlabel('Time (days)');
% hold on
% ylim([-1000 800]);
% ax.FontSize = p.font_size;
% yyaxis right
% plot(t/86400, u.T_CD(t),'color',[0.6350 0.0780 0.1840], 'LineWidth',1); 
% ylabel('T_C_D (^oC)');
% ylim([2 32]);
% hold on
% ax.FontSize = p.font_size;
% hold off

%% Define Process Parameters
% General
p.rho_Water = 1000;      % kg/m3,  Density of water
p.h_0       = 0.10186;   % kJ/kg,  Reference specific enthalpy
p.T_0       = 0.01;      % oC,     Reference temperature
p.C_p       = 4.1831;    % kJ/kgC, Heat capacity of water

% Square Dam
p.m_SDmax   = 6000000;   % kg,   Maximum mass capacity of SD
p.height_SD = 3;         % m,    Height of the SD
p.area_SD   = 2000;      % m2,   Surface area of the SD
p.m_evapSD =  1;         % kg/s, Rate of evaporation from SD surface

% Pre-cooling Towers
p.rho_Air   = 1.225;     % kg/m3,  Density of air
p.H_outAir  = 100;       % %,      Humidity of air leaving towers
p.A_enthal  = 1.01;      % -,      Enthalpy equation constant
p.B_enthal  = 1.89;      % -,      Enthalpy equation constant
p.C_enthal  = 2501;      % -,      Enthalpy equation constant
p.C1        = -5.8E+3;   % -,      Pressure equation constant
p.C2        = 1.391;     % -,      Pressure equation constant
p.C3        = -4.864E-2; % -,      Pressure equation constant
p.C4        = 4.176E-5;  % -,      Pressure equation constant
p.C5        = -1.445E-8; % -,      Pressure equation constant
p.C6        = 6.546;     % -,      Pressure equation constant
p.M_r       = 0.62198;   % -,      Ratio between molar mass water and air
p.P_Tot     = 101.325;   % kPa,    Atmospheric pressure
p.height_PT = 3;         % m, Height of the Pre-cooling Tower basins
p.length_PT = 14.8;      % m, Length of the Pre-cooling Tower basins
p.width_PT  = 10.45*2;   % m, Width of the Pre-cooling Tower basins
p.area_PT   = p.length_PT*p.width_PT; % m2, Area of the Pre-cooling Tower basins
p.volume_PT = p.area_PT*p.height_PT;  % m3, Volume of the Pre-cooling Tower basins
p.m_PTAmax  = p.volume_PT/2*1000;     % kg, Maximum mass capacity of PT A basin
p.m_PTBmax  = p.volume_PT/2*1000;     % kg, Maximum mass capacity of PT B basin
p.m_PTmax   = p.m_PTAmax*2;
p.m_Air     = 800;                    % kg/s, Estimated mass flowrate of air

% Fridge Plants
p.m_RPj   = 10000;        % kg,     Mass held in the evaporator of each fridge plant
p.UA_RP1  = 15000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 1
p.UA_RP2  = 15000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 2
p.UA_RP3  = 15000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 3
p.UA_RP4  = 32500;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 4
p.UA_RP5  = 23000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 5
p.UA_RP   = [p.UA_RP1, p.UA_RP2, p.UA_RP3, p.UA_RP4, p.UA_RP5];
p.UA_amb1 = -10;          % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 1
p.UA_amb2 = -10;          % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 2
p.UA_amb3 = -300;          % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 3
p.UA_amb4 = -300;          % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 4
p.UA_amb5 = -50;          % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 5
p.UA_amb  = [p.UA_amb1, p.UA_amb2, p.UA_amb3, p.UA_amb4, p.UA_amb5];

% Chill Dam
p.m_CDmax   = 8000000; % kg, Maximum mass capacity of CD
p.height_CD = 6.7;     % m,  Height of the Chill Dam
p.length_CD = 48;      % m,  Length of the Chill Dam
p.width_CD  = 25;      % m,  Width of the Chill Dam
p.area_CD   = p.length_CD*p.width_CD; % m2, Area of the Chill Dam
p.UA_CD     = 190;     % kJ/Ks, Heat transfer coefficient for the rate
                       % of heat transfer between the water in the CD and
                       % the environment through the dam walls

% Electricity Tariffs
p.offpeak_tariff  = 0.45; % R/kWh, Minimum or off-peak hourly charge for electricity by Eskom
p.standard_tariff = 0.83; % R/kWh, Intermediate or standard hourly charge for electricity by Eskom
p.onpeak_tariff   = 2.76; % R/kWh, Maximum or on-peak hourly charge for electricity by Eskom

% Electrical Tariff Vectors
p.weekday_tariff_vec  = [ones(1,6)*p.offpeak_tariff ones(1,3)*p.onpeak_tariff ones(1,8)*p.standard_tariff...
                        ones(1,2)*p.onpeak_tariff ones(1,3)*p.standard_tariff ones(1,2)*p.offpeak_tariff];
p.saturday_tariff_vec = [ones(1,7)*p.offpeak_tariff ones(1,5)*p.standard_tariff...
                        ones(1,6)*p.offpeak_tariff ones(1,2)*p.standard_tariff ones(1,4)*p.offpeak_tariff];
p.sunday_tariff_vec   = [ones(1,24)*p.offpeak_tariff];

p.weekday_tariff_vec_min  = [ones(1,6*60)*p.offpeak_tariff ones(1,3*60)*p.onpeak_tariff ones(1,8*60)*p.standard_tariff...
                            ones(1,2*60)*p.onpeak_tariff ones(1,3*60)*p.standard_tariff ones(1,2*60)*p.offpeak_tariff];
p.saturday_tariff_vec_min = [ones(1,7*60)*p.offpeak_tariff ones(1,5*60)*p.standard_tariff...
                            ones(1,6*60)*p.offpeak_tariff ones(1,2*60)*p.standard_tariff ones(1,4*60)*p.offpeak_tariff];
p.sunday_tariff_vec_min   = [ones(1,24*60)*p.offpeak_tariff];

p.whole_week_tariff_vec_minutes = [repmat(p.weekday_tariff_vec_min,1,5) p.saturday_tariff_vec_min p.sunday_tariff_vec_min];

% Power Pull by Fridge Plants
p.pull_full_efficiency = 3000; % kW, Power pulled by ONE fridge plant when operating at 100% efficiency
p.onpeak_consumption  = 3500;  % kW, Maximum power consumption allowed during peak electricity times
p.offpeak_consumption = 16000; % kW, Maximum power consumption allowed during off peak electricity times

% Power Pull Limit Vectors
p.weekday_power_vec  = [ones(1,6)*p.offpeak_consumption ones(1,3)*p.onpeak_consumption ones(1,8)*p.offpeak_consumption...
                        ones(1,2)*p.onpeak_consumption ones(1,3)*p.offpeak_consumption ones(1,2)*p.offpeak_consumption];
p.saturday_power_vec = [ones(1,24)*p.offpeak_consumption];
p.sunday_power_vec   = [ones(1,24)*p.offpeak_consumption];
%% Load AR Model Data
% Load the AR model data, providing stochastic sets of data
% representing the trends in the measurement data, with
% calculated variances and constants.
[p,u] = ARModelsSD(t,u,p);
[p,u] = ARModelsPT(t,u,p); 
[p,u] = ARModelsRP(t,u,p);
[p,u] = ARModelsCD(t,u,p);

%% Define Kalman Filter parameters
% Square Dam
p.xhat_0_SD = [2.9; 400];               % Intial state estimates (Level in m, inlet flowrate in L/s)
p.H_SD      = eye(length(p.xhat_0_SD)); % Observation matrix 
p.c_K_SD = [0 ; p.c_F_in_SD];           % Constants from AR model
p.w_F_in_SD = 0.1;                      % Tuning for F_in_SD process noise
p.w_K_SD = [p.w_L_SD ; p.w_F_in_SD];    % Variance from AR model
p.L_sensor_noise_SD   = 0.01;           % Std dev of measurement noise, level sensor error
p.F_in_meter_noise_SD = 0.01;           % Std dev of measurement noise, inlet stream flowmeter error
p.R_SD = diag([p.L_sensor_noise_SD^2 ;...
            p.F_in_meter_noise_SD^2]);  % Measurement noise covariance matrix (measurement error squared)   
p.Q_SD = diag([0.09 ; p.w_K_SD(2)]);    % Process noise covariance matrix based on variance from AR model 
p.A_SD = [0 p.C_L_SD;    
       0 0];                            % Transition matrix

% Pre-cooling Towers
p.xhat_0_PT = [2.9; 10; 15; 25; 13; 20; 50]; % Intial state estimates (Level in m, inlet flowrate in L/s)
p.H_PT      = eye(length(p.xhat_0_PT)); % Observation matrix 
p.c_K_PT = [0 ; p.c_F_Antisurge; p.c_T_Antisurge; p.c_T_out_SD; 0; 0; p.c_H_Env]; % Constants from AR model
p.w_F_Antisurge = 0.000000000001;          % Tuning for F_Antisurge process noise
p.w_T_PT = 0.0001;                      % Tuning for T_PT process noise
p.w_H_Env = 0.1;                        % Tuning for H_Env process noise
p.w_K_PT = [p.w_L_PT ; p.w_F_Antisurge; p.w_T_Antisurge; p.w_T_out_SD; p.w_T_PT; 0.002; p.w_H_Env];    % Variance from AR model
p.L_sensor_noise_PT   = 0.01;           % Std dev of measurement noise, level sensor error
p.F_meter_noise_PT = 0.01;           % Std dev of measurement noise, inlet stream flowmeter error
p.T_sensor_noise_PT = 0.01;           % Std dev of measurement noise, inlet stream flowmeter error
p.H_sensor_noise_PT = 0.01;           % Humidity sensor noise
p.R_PT = diag([p.L_sensor_noise_PT^2 ;...
            p.F_meter_noise_PT^2;...
            p.T_sensor_noise_PT^2 ;...
            p.T_sensor_noise_PT^2 ;...
            p.T_sensor_noise_PT^2 ;...
            p.T_sensor_noise_PT^2 ;...
            p.H_sensor_noise_PT^2]);     % Measurement noise covariance matrix (measurement error squared)  
p.Q_PT = diag([p.w_K_PT(1); p.w_K_PT(2); p.w_K_PT(3); p.w_K_PT(4);...
               p.w_K_PT(5); p.w_K_PT(6); p.w_K_PT(7)]); % Process noise covariance matrix based on variance from AR model 


% Fridge Plants
p.xhat_0_RP = [4; 1; -0.05; 12];        % Intial state estimates: T_RPj, Q_refrj, Q_amb, T_inRPj
p.H_RP      = eye(length(p.xhat_0_RP)); % Observation matrix 
p.c_K_RP = [0 0 0 0 0; p.c_Q_refr; p.c_Q_amb; p.c_T_in_RP]; % Constants from AR model
%p.w_T_in_RP = [0.001 0.001 0.004 0.003 0.002]; % Tuning for T_in_RP process noise
p.w_K_RP = [p.w_T_RP ; p.w_Q_refr; p.w_Q_amb; p.w_T_in_RP]; % Variance from AR model
p.T_sensor_noise_RP   = 0.01;         % Std dev of measurement noise, level sensor error
p.Q_lumped_noise      = 0.01;         % Std dev of measurement noise, inlet stream flowmeter error
p.R_RP = diag([p.T_sensor_noise_RP^2 ;...
            p.Q_lumped_noise^2;...
            p.Q_lumped_noise^2 ;...
            p.T_sensor_noise_RP^2 ]); % Measurement noise covariance matrix (measurement error squared) 
p.Q1 = diag([p.w_T_RP(1); 0.002;  1E-12; 1E-12]); % Process noise covariance matrix based on variance from AR model for RP1
p.Q2 = diag([p.w_T_RP(2); 0.02;  1E-12; 1E-12]);  % Process noise covariance matrix based on variance from AR model for RP2
p.Q3 = diag([p.w_T_RP(3); 0.02;  1E-12; 1E-12]);  % Process noise covariance matrix based on variance from AR model for RP3
p.Q4 = diag([p.w_T_RP(4); 0.02;  1E-12; 1E-12]);  % Process noise covariance matrix based on variance from AR model for RP4
p.Q5 = diag([p.w_T_RP(5); 0.02;  1E-12; 1E-12]);  % Process noise covariance matrix based on variance from AR model for RP5

% Chill Dam
p.xhat_0_CD = [2.9; 4; 300; 300];        % Intial state estimates (Level in m, inlet flowrate in L/s)
p.H_CD      = eye(length(p.xhat_0_CD));  % Observation matrix 
p.c_K_CD = [0 ; 0; p.c_F_Ice; p.c_F_UG]; % Constants from AR model
p.w_F_Ice = 0.00000000000000000001;
p.w_F_UG = 0.001;
p.w_T_CD = 0.00001;                      % Tuning for T_CD process noise
p.w_K_CD = [p.w_L_CD ; p.w_T_CD ; p.w_F_Ice ; p.w_F_UG]; % Variance from AR model
p.L_sensor_noise_CD   = 0.01;           % Std dev of measurement noise, level sensor error
p.T_sensor_noise_CD   = 0.01;
p.F_meter_noise_CD    = 0.01;           % Std dev of measurement noise, inlet stream flowmeter error
p.R_CD = diag([p.L_sensor_noise_CD^2 ;...
               p.T_sensor_noise_CD^2 ;...
               p.F_meter_noise_CD^2  ;...
               p.F_meter_noise_CD^2]);  % Measurement noise covariance matrix (measurement error squared)
p.Q_CD = diag([0.009; p.w_K_CD(2); p.w_K_CD(3) ; p.w_K_CD(4)]); % Process noise covariance matrix based on variance from AR model 

%% Define MPC parameters
% General
p.N          = 3;                           % ~, Number of samples of the control input/prediction nodes
p.N_Econ     = 1;                           % ~, Prediction horizon for optimiser/how often control actions can take place
p.Ts         = 60;                          % s, Sampling period, the frequency at which a new control input is determined
p.Stp        = length(t);                   % ~, Number of steps in simulation
p.TL         = (p.Ts*p.Stp) - p.Ts;         % s, Total time or Time Limit (sampling period x total time)
%p.loop       = 2:1:720;                    % ~, Simulation points for the loop in minutes (mostly for testing)
p.loop       = 2:1:2879                    % 48 hours
p.loop_hours = 1:1:24;                      % ~, Simulation points for the loop in hours (mostly for testing)
p.SP_changes = 2296;                        % ~, Number of SP changes (2296 - every 22 min. 1335 - every 38 min)
p.SP_times   = (0:p.TL/p.SP_changes:p.TL)'; % Times at which the SP should change
p.font_size  = 18;

% Square Dam
p.L_SS_SD       = 70;                          % %, Steady state level in the square dam (initial PV value)
p.F_outSS_SD    = u.F_outSD(1);                % L/s, Steady state outlet flowrate (initial MV value)
p.uvec_init_SD  = p.F_outSS_SD*ones(1,p.N);    % L/s, Initial points (sequence guess)
p.SP_SD         = 85*ones(1,p.N);              % %, Initial SP for the level in the square dam (initial SP value)
p.SP_min_SD     = 75;                          % %, Lowest SP for the level in the Dam
p.SP_max_SD     = 95;                          % %, Highest SP for the level in the Dam
p.SP_samples_SD = p.SP_min_SD...
                + (p.SP_max_SD - p.SP_min_SD)...
                * rand(p.SP_changes, 1);       % Sample SP changes
p.MV_min_SD     = 0*ones(1,p.N);               % L/s, Minimum MV limit
p.MV_max_SD     = 750*ones(1,p.N);             % L/s, Maximum MV limit
p.PV_min_SD     = 40*ones(1,p.N);              % %, Minimum PV limit (minimum level)
p.PV_max_SD     = 95*ones(1,p.N);              % %, Maximum PV limit (maximum level)
p.Q_Weight_SD   = 0.1;                         % SP weight
p.R_Weight_SD   = 1;                           % MV weight

% Pre-cooling Towers
p.L_SS_PT       = 70;                          % %, Steady state level in the pre-cooling towers (initial PV value
p.F_outSS_PT    = u.F_outPT(1);                % L/s, Steady state outlet flowrate (initial MV value)
p.uvec_init_PT  = p.F_outSS_PT*ones(1,p.N);    % L/s, Initial points (sequence guess)
p.SP_PT         = 85*ones(1,p.N);              % %, Initial SP for the level in the square dam (initial SP value)
p.SP_min_PT     = 75;                          % %, Lowest SP for the level in the basins
p.SP_max_PT     = 90;                          % %, Highest SP for the level in the basins
p.SP_samples_PT = p.SP_min_PT...
                + (p.SP_max_PT - p.SP_min_PT)...
                * rand(p.SP_changes, 1);       % Sample SP changes
p.MV_min_PT     = 250*ones(1,p.N);             % L/s, Minimum MV limit (ensure that there is always enough flow for RPs)
p.MV_max_PT     = 700*ones(1,p.N);             % L/s, Maximum MV limit
p.PV_min_PT     = 40*ones(1,p.N);              % %, Minimum PV limit (minimum level)
p.PV_max_PT     = 98*ones(1,p.N);              % %, Maximum PV limit (maximum level)
p.Q_Weight_PT   = 1;                           % SP weight
p.R_Weight_PT   = 0.1;                         % MV weight

% Fridge Plants
p.T_SS_RP       = 13;                          % oC, Steady state inlet temperature of the fridge plants (initial PV value)
p.F_RecSS_RP    = 50;                          % L/s, Steady state recycle flowrate (initial MV value)
p.uvec_init_RP  = p.F_RecSS_RP*ones(1,p.N);    % L/s, Initial points (sequence guess)
p.SP_RP         = 11*ones(1,p.N);              % oC, Initial SP for the outlet temperature of the fridge plants (initial SP value)
p.SP_min_RP     = 10;                          % oC, Lowest SP for the outlet temperature of the fridge plants
p.SP_max_RP     = 12;                          % oC, Highest SP for the outlet temperature of the fridge plants
p.SP_samples_RP = p.SP_min_RP...
                + (p.SP_max_RP - p.SP_min_RP)...
                * rand(p.SP_changes, 1);       % Sample SP changes
p.MV_min_RP     = 0*ones(1,p.N);               % L/s, Minimum MV limit (ensure that there is always enough flow for RPs)
p.MV_max_RP     = 150*ones(1,p.N);             % L/s, Maximum MV limit
p.PV_min_RP     = 10*ones(1,p.N);              % oC, Minimum PV limit (minimum temperature)
p.PV_max_RP     = 12*ones(1,p.N);              % oC, Maximum PV limit (maximum temperature)
p.Q_Weight_RP   = 5;                           % SP weight
p.R_Weight_RP   = 0.1;                         % MV weight

% Chill Dam
p.L_SS_CD       = 70;                          % %, Steady state level in the chill dam (initial PV value)
p.T_in_SS_RP    = 13;                          % oC, Steady state temperature entering the RPs (initial PV value)
p.F_outSS_CD    = u.F_outCD(1);                % L/s, Steady state outlet flowrate (initial MV value)
p.uvec_init_CD  = p.F_outSS_CD*ones(1,p.N);    % L/s, Initial points (sequence guess)
p.SP_L_CD       = 85*ones(1,p.N);              % %, Initial SP for the level in the .chill dam (initial SP value)
p.SP_T_CD       = 12*ones(1,p.N);              % %, Initial SP for the temperature entering the RPs (initial SP value)
p.SP_L_min_CD   = 40;                          % %, Lowest SP for the level in the Dam
p.SP_L_max_CD   = 95;                          % %, Highest SP for the level in the Dam
p.SP_T_min_CD   = 10;                          % %, Lowest SP for the temperature entering the RPs
p.SP_T_max_CD   = 13;                          % %, Highest SP for the temperature entering the RPs
p.SP_samples_L_CD = p.SP_L_min_CD...
                + (p.SP_L_max_CD - p.SP_L_min_CD)...
                * rand(p.SP_changes, 1);       % Sample SP changes for level
p.SP_samples_T_CD = p.SP_T_min_CD...
                + (p.SP_T_max_CD - p.SP_T_min_CD)...
                * rand(p.SP_changes, 1);       % Sample SP changes for temperature
p.MV_min_CD     = 0*ones(1,p.N);               % L/s, Minimum MV limit (ensure that there is always enough flow for RPs)
p.MV_max_CD     = 750*ones(1,p.N);             % L/s, Maximum MV limit
p.PV_min_L_CD   = 15*ones(1,p.N);              % %, Minimum PV limit (minimum level)
p.PV_max_L_CD   = 95*ones(1,p.N);              % %, Maximum PV limit (maximum level)
p.PV_min_T_CD   = 10*ones(1,p.N);              % %, Minimum PV limit (minimum level)
p.PV_max_T_CD   = 12.5*ones(1,p.N);              % %, Maximum PV limit (maximum level)
p.Q_Weight_CD   = 1;                           % SP weight
p.R_Weight_CD   = 0.1;                         % MV weight

% Economic Optimisation
p.Econ_SS        = 3000;                       % kW, Steady state energy consumption (initial PV value)
p.L_SD_SP_SS     = 80;                         % %, Steady state Square Dam Level SP (initial MV value)
p.L_PT_SP_SS     = 80;                         % %, Steady state Pre-Cooling Tower Level SP (initial MV value)
p.L_CD_SP_SS     = 80;                         % %, Steady state Chill Dam Level SP (initial MV value)
p.T_RP_SP_SS     = 11;                         % oC, Steady state Fridge Plant Inlet Temp SP (initial MV value)
p.T_CD_SP_SS     = 12;                         % oC, Steady state Chill Dam Outlet Temp SP (initial MV value)
p.no_RP_SP_SS    = 3;                          % #, Steady state Number of RPs operating SP (initial MV value)
p.uvec_init_Econ = p.no_RP_SP_SS*ones(1,p.N_Econ);  % -, Initial points for MV (sequence guess)
p.MV_min_Econ    = 0*ones(1,p.N_Econ);              % -, Minimum MVs limit for number of RPs operating
p.MV_max_Econ    = 5*ones(1,p.N_Econ);              % L/s, Maximum MVs limit for number of RPs operating
p.Q_Weight_Econ   = 1;                         % SP weight
p.R_Weight_Econ   = 0.1;                       % MV weight
%% Define state structure and initial conditions

% Square Dam
s.statefields_SD    = {'L_SD'}; 
s.MPCstatefields_SD = {'L_SD','F_in_SD'};
s.KFstatefields_SD  = {'L_SD', 'F_in_SD', 'P_SD'}; 
x0.L_SD    = u.L_SD(0);                      
x0.F_in_SD = u.F_inSD(0);                
x0.P_SD    = 0.1;

% Pre-cooling Towers
s.statefields_PT    = {'L_PT','T_PT'};                
s.MPCstatefields_PT = {'L_PT','F_AntiSurge','T_AntiSurge','T_outSD','T_PT', 'T_amb','H_Env'};      
s.KFstatefields_PT  = {'L_PT','F_AntiSurge','T_AntiSurge','T_outSD','T_PT', 'T_amb','H_Env','P_PT'}; 
x0.L_PT    = (u.L_PTA(0) + u.L_PTB(0))/2;     
x0.F_AntiSurge = u.F_AntiSurge(0);
x0.T_AntiSurge = u.T_AntiSurge(0);
x0.T_outSD = u.T_outSD(0);
x0.T_PT  = 15;  
x0.T_amb = 20;
x0.H_Env = 50;
x0.P_PT  = 0.1; 


% Fridge Plants
s.statefields_RP    = {'T1',      'T2',      'T3',      'T4',      'T5'};
s.MPCstatefields_RP = {'T1',      'T2',      'T3',      'T4',      'T5';...
                       'Qrefr1',  'Qrefr2',  'Qrefr3',  'Qrefr4',  'Qrefr5';...
                       'Qamb1',   'Qamb2',   'Qamb3',   'Qamb4',   'Qamb5';...
                       'Tin1',    'Tin2',    'Tin3',    'Tin4',    'Tin5'};
s.KFstatefields_RP  = {'T1',      'T2',      'T3',      'T4',      'T5';...
                       'Qrefr1',  'Qrefr2',  'Qrefr3',  'Qrefr4',  'Qrefr5';...
                       'Qamb1',   'Qamb2',   'Qamb3',   'Qamb4',   'Qamb5';...
                       'Tin1',    'Tin2',    'Tin3',    'Tin4',    'Tin5';...
                       'P1',      'P2',      'P3',      'P4',      'P5'};   
x0.T1 = u.T_RP.Values(1,1);  
x0.T2 = u.T_RP.Values(1,2);
x0.T3 = u.T_RP.Values(1,3);
x0.T4 = u.T_RP.Values(1,4);
x0.T5 = u.T_RP.Values(1,5);
x0.Qrefr1 = u.Q_refr_generated.Values(1,1);  
x0.Qrefr2 = u.Q_refr_generated.Values(1,2);
x0.Qrefr3 = u.Q_refr_generated.Values(1,3);
x0.Qrefr4 = u.Q_refr_generated.Values(1,4);
x0.Qrefr5 = u.Q_refr_generated.Values(1,5);
x0.Qamb1 = u.Q_amb_generated.Values(1,1);  
x0.Qamb2 = u.Q_amb_generated.Values(1,2);
x0.Qamb3 = u.Q_amb_generated.Values(1,3);
x0.Qamb4 = u.Q_amb_generated.Values(1,4);
x0.Qamb5 = u.Q_amb_generated.Values(1,5);
x0.Tin1 = u.T_in_generatedRP.Values(1,1);  
x0.Tin2 = u.T_in_generatedRP.Values(1,2);
x0.Tin3 = u.T_in_generatedRP.Values(1,3);
x0.Tin4 = u.T_in_generatedRP.Values(1,4);
x0.Tin5 = u.T_in_generatedRP.Values(1,5);
x0.P_RP = 0.1; 

% Chill Dam
s.statefields_CD    = {'L_CD','T_CD'};              
s.MPCstatefields_CD = {'L_CD','T_CD','F_Ice','F_UG'};    
s.KFstatefields_CD  = {'L_CD','T_CD','F_Ice','F_UG','P_CD'};
x0.L_CD = u.L_CD(0);   
x0.T_CD = 4;         
x0.F_Ice = 100; 
x0.F_UG  = 200;
x0.P_CD = 0.1;        

save CombinedSystem.mat s x0 p u t
%% Breakpoint for testing
clc
clear
load CombinedSystem.mat
%% MPC Initialisation

% General
options = optimoptions('fmincon','Display','off');
response.CD = [x0.L_CD; x0.T_CD];
response.RP = [4 4 4 4 4];
output.MV_CD = @(t) p.F_outSS_CD;
output.MV1   = @(t) p.F_RecSS_RP;
output.MV2   = @(t) p.F_RecSS_RP;
output.MV3   = @(t) p.F_RecSS_RP;
output.MV4   = @(t) p.F_RecSS_RP;
output.MV5   = @(t) p.F_RecSS_RP;
output.MV_L_SD_SP = @(t) p.L_SD_SP_SS;
output.MV_L_PT_SP = @(t) p.L_PT_SP_SS;
output.MV_L_CD_SP = @(t) p.L_CD_SP_SS;
output.MV_T_RP_SP = @(t) p.T_RP_SP_SS;
output.MV_T_CD_SP = @(t) p.T_CD_SP_SS;
output.MV_no_RP_SP = @(t) p.no_RP_SP_SS;
output.MV_Econ = ones(p.N_Econ,1).*output.MV_no_RP_SP(t);

% Square Dam
solSD.y = p.L_SS_SD;        
z.L_SD = []; z.F_in_SD = [];


% Pre-cooling Towers
solPT.y = [p.L_SS_PT; x0.T_PT]; 
z.L_PT = [];    z.F_AntiSurge = []; z.T_AntiSurge = []; 
z.T_outSD = []; z.T_PT = []; z.T_amb = []; z.H_Env = [];


% Fridge Plants
sol1.y = [x0.T1];
sol2.y = [x0.T2];
sol3.y = [x0.T3];
sol4.y = [x0.T4];
sol5.y = [x0.T5];
z.T1     = [];     z.T2 = [];     z.T3 = [];     z.T4 = [];     z.T5 = [];
z.Qrefr1 = []; z.Qrefr2 = []; z.Qrefr3 = []; z.Qrefr4 = []; z.Qrefr5 = []; 
z.Qamb1  = [];  z.Qamb2 = [];  z.Qamb3 = [];  z.Qamb4 = [];  z.Qamb5 = []; 
z.Tin1   = [];   z.Tin2 = [];   z.Tin3 = [];   z.Tin4 = [];   z.Tin5 = [];
v.T_inRP = [p.T_SS_RP p.T_SS_RP p.T_SS_RP p.T_SS_RP p.T_SS_RP];


% Chill Dam
solCD.y = [p.L_SS_CD; x0.T_CD]; 
z.L_CD = []; z.T_CD = []; z.F_Ice = []; z.F_UG = [];


% Generate measurements for all units
z = Meas(solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, u, 0, p, z, v);   
x = z;   

% Add to the x structure
x.P_SD = x0.P_SD;
x.P_PT = 0.1;
x.P_CD = 0.1;
x.P1 = 0.1; x.P2 = 0.1; x.P3 = 0.1; x.P4 = 0.1; x.P5 = 0.1;

% Square Dam Optimisation
u_opt_SD = fmincon(@(uMV) costSD(t(1), uMV, u, p, s, x, output), p.uvec_init_SD,...
    [], [], [], [], p.MV_min_SD, p.MV_max_SD, [], options);                                                     
MV_SD = u_opt_SD(1);
output.MV_SD = @(t) MV_SD; 
solSD = ode45(@(t,x) SquareDamODEs(s,p,x,u,t,output,z), [t(1) t(2)], solSD.y); 
response.SD(1,:) = deval(solSD, t(1));

% Pre-cooling Towers Optimisation
u_opt_PT = fmincon(@(uMV) costPT(t(1), uMV, u, p, s, x, output), p.uvec_init_PT,...
    [], [], [], [], p.MV_min_PT, p.MV_max_PT, [], options);                                                   
MV_PT = u_opt_PT(1); 
output.MV_PT = @(t) MV_PT; 
solPT = ode45(@(t,x) PreCoolingTowersODEs(s,p,x,u,t,output,z), [t(1) t(2)], solPT.y); 
response.PT = deval(solPT(1,:), t(1));


% Fridge Plants Optimisation
% Fridge Plant 1
u_opt1 = fmincon(@(uMV) cost1(t(1), uMV, u, p, s, x, 1, output, response), p.uvec_init_RP,...
         [], [], [], [],p.MV_min_RP, p.MV_max_RP, [], options); 
MV1 = u_opt1(1);
output.MV1 = @(t) MV1; 
sol1 = ode45(@(t,x) FridgePlantsODEs1(s,p,x,u,t,output,1,response), [t(1) t(2)], sol1.y); 
response.RP1(1,:) = deval(sol1, t(1));

% Fridge Plant 2
u_opt2 = fmincon(@(uMV) cost2(t(1), uMV, u, p, s, x, 2, output, response), p.uvec_init_RP,...
         [], [], [], [], p.MV_min_RP, p.MV_max_RP, [], options); 
MV2 = u_opt2(1);
output.MV2 = @(t) MV2; 
sol2 = ode45(@(t,x) FridgePlantsODEs2(s,p,x,u,t,output,2,response), [t(1) t(2)], sol2.y); 
response.RP2(1,:) = deval(sol2, t(1));

% Fridge Plant 3
u_opt3 = fmincon(@(uMV) cost3(t(1), uMV, u, p, s, x, 3, output, response), p.uvec_init_RP,...
         [], [], [], [], p.MV_min_RP, p.MV_max_RP, [], options); 
MV3 = u_opt3(1);
output.MV3 = @(t) MV3; 
sol3 = ode45(@(t,x) FridgePlantsODEs3(s,p,x,u,t,output,3,response), [t(1) t(2)], sol3.y); 
response.RP3(1,:) = deval(sol3, t(1));

% Fridge Plant 4
u_opt4 = fmincon(@(uMV) cost4(t(1), uMV, u, p, s, x, 4, output, response), p.uvec_init_RP,...
         [], [], [], [], p.MV_min_RP, p.MV_max_RP, [], options); 
MV4 = u_opt4(1);
output.MV4 = @(t) MV4; 
sol4 = ode45(@(t,x) FridgePlantsODEs4(s,p,x,u,t,output,4,response), [t(1) t(2)], sol4.y);
response.RP4(1,:) = deval(sol4, t(1));

% Fridge Plant 5
u_opt5 = fmincon(@(uMV) cost5(t(1), uMV, u, p, s, x, 5, output, response), p.uvec_init_RP,...
         [], [], [], [], p.MV_min_RP, p.MV_max_RP, [], options); 
MV5 = u_opt5(1);
output.MV5 = @(t) MV5; 
sol5 = ode45(@(t,x) FridgePlantsODEs5(s,p,x,u,t,output,5,response), [t(1) t(2)], sol5.y); 
response.RP5(1,:) = deval(sol5, t(1));

response.RP = [response.RP1 response.RP2 response.RP3 response.RP4 response.RP5];


% Chill Dam Optimisation
u_opt_CD = fmincon(@(uMV) costCD(t(1), uMV, u, p, s, x, response, output), p.uvec_init_CD,...
    [], [], [], [], p.MV_min_CD, p.MV_max_CD, [], options);                                                 
MV_CD = u_opt_CD(1); 
output.MV_CD = @(t) MV_CD; 
solCD = ode45(@(t,x) ChillDamODEs(s,p,x,u,t,output,response), [t(1) t(2)], solCD.y); % Calculate the ground truth for the first time interval
response.CD = deval(solCD(1,:), t(1));

% Scheduling Optimisation (Economic Optimisation)
% GA_options = optimoptions('ga','Display','off');
% u_opt_Econ = fmincon(@(uMV) costEcon(t(1), uMV, u, p, s, x, output), p.uvec_init_Econ,...
%     [], [], [], [], p.MV_min_Econ, p.MV_max_Econ, [], options);   
MV_Econ = output.MV_Econ;
% output.MV_Econ = MV_Econ; 

            % GA_options = optimoptions('ga','Display','off', 'MaxGenerations', 1);
            % u_opt_Econ = ga(@(uMV) costEcon(t(2), uMV, u, p, s, x, output), 1,...
            %      [], [], [], [], p.MV_min_Econ, p.MV_max_Econ, ...
            %      @(uMV) mycon(t(2), uMV, x, output, s, u, z, p, response, t, ...
            %      solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, ...
            %      u_opt_SD, u_opt_PT, u_opt1, u_opt2, u_opt3, u_opt4, u_opt5, u_opt_CD),...
            %      1, GA_options);   
            % output.MV_Econ(end+1) = u_opt_Econ(1);

% Save the generated set-points
% for b = 1:2
%     saved.SP_L_SD(b,:) = output.MV_Econ(1,end); 
%     saved.SP_L_PT(b,:) = output.MV_Econ(2,end);
%     saved.SP_L_CD(b,:) = output.MV_Econ(3,end);
%     saved.SP_T_RP(b,:) = output.MV_Econ(4,end);   
%     saved.SP_T_CD(b,:) = output.MV_Econ(5,end);
%     saved.SP_no_RP(b,:) = output.MV_Econ(6,end);
% end
for b = 1:2
    saved.SP_L_SD(b,:) = p.L_SD_SP_SS; 
    saved.SP_L_PT(b,:) = p.L_PT_SP_SS;
    saved.SP_L_CD(b,:) = p.L_CD_SP_SS;
    saved.SP_T_RP(b,:) = p.T_RP_SP_SS;   
    saved.SP_T_CD(b,:) = p.T_CD_SP_SS;
    saved.SP_no_RP(b,:) = output.MV_Econ(end);
end

T_UG = 3;

%% MPC Loop

    for j = p.loop
        
    v = RPIntermediatesKF(x, u, p, t(j), output, s, response);
    z = Meas(solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, u, t(j), p, z, v);
    x = KalmanFilter(x, s, output, u, z, [t(j-1) t(j)], p, response);


    % Square Dam
    u_opt_SD = fmincon(@(uMV) costSD(t(j), uMV, u, p, s, x, output), u_opt_SD, [], [], [], [],...
               p.MV_min_SD, p.MV_max_SD, [], options);
    MV_SD(end+1) = u_opt_SD(1); 
    output.MV_SD = griddedInterpolant(t(1:j), MV_SD, 'previous');
    solSD = odextend(solSD, @(t,x) SquareDamODEs(s, p, x, u, t, output,z), t(j+1)); 
    response.SD(j,:) = deval(solSD, t(j));

    % Pre-cooling Towers
    p.MV_min_PT = 50 + output.MV_Econ(end)*50*ones(1,p.N); % Determine the lower MV limits based on number of fridge plants operating
    % p.MV_max_PT = ((output.MV_Econ(end)*200)+200)*ones(1,p.N); % Determine the upper MV limits based on number of fridge plants operating
    u_opt_PT = fmincon(@(uMV) costPT(t(j), uMV, u, p, s, x, output), u_opt_PT, [], [], [], [],...
               p.MV_min_PT, p.MV_max_PT, [], options);
    MV_PT(end+1) = u_opt_PT(1); 
    output.MV_PT = griddedInterpolant(t(1:j), MV_PT, 'previous');
    solPT = odextend(solPT, @(t,x) PreCoolingTowersODEs(s, p, x, u, t, output,z), t(j+1));
    response.PT(:,j) = deval(solPT, t(j));
    
    % Fridge Plant 1
    u_opt1 = fmincon(@(uMV1) cost1(t(j), uMV1, u, p, s, x, 1, output, response), u_opt1,...
         [], [], [], [], p.MV_min_RP, p.MV_max_RP, [], options); 
    MV1(end+1) = u_opt1(1);
    output.MV1 = griddedInterpolant(t(1:j), MV1, 'previous');
    sol1 = odextend(sol1, @(t,x) FridgePlantsODEs1(s, p, x, u, t, output, 1,response), t(j+1)); 
    response.RP1(j,:) = deval(sol1, t(j));
   
    % Fridge Plant 2
    u_opt2 = fmincon(@(uMV2) cost2(t(j), uMV2, u, p, s, x, 2, output, response), u_opt2,...
         [], [], [], [], p.MV_min_RP, p.MV_max_RP, [], options); 
    MV2(end+1) = u_opt2(1);
    output.MV2 = griddedInterpolant(t(1:j), MV2, 'previous');
    sol2 = odextend(sol2, @(t,x) FridgePlantsODEs2(s, p, x, u, t, output, 2,response), t(j+1)); 
    response.RP2(j,:) = deval(sol2, t(j));

    % Fridge Plant 3
    u_opt3 = fmincon(@(uMV3) cost3(t(j), uMV3, u, p, s, x, 3, output, response), u_opt3,...
         [], [], [], [], p.MV_min_RP, p.MV_max_RP, [], options); 
    MV3(end+1) = u_opt3(1);
    output.MV3 = griddedInterpolant(t(1:j), MV3, 'previous');
    sol3 = odextend(sol3, @(t,x) FridgePlantsODEs3(s, p, x, u, t, output, 3,response), t(j+1));
    response.RP3(j,:) = deval(sol3, t(j));

    % Fridge Plant 4
    u_opt4 = fmincon(@(uMV4) cost4(t(j), uMV4, u, p, s, x, 4, output, response), u_opt4,...
         [], [], [], [], p.MV_min_RP, p.MV_max_RP, [], options); 
    MV4(end+1) = u_opt4(1);
    output.MV4 = griddedInterpolant(t(1:j), MV4, 'previous');
    sol4 = odextend(sol4, @(t,x) FridgePlantsODEs4(s, p, x, u, t, output, 4,response), t(j+1));
    response.RP4(j,:) = deval(sol4, t(j));

    % Fridge Plant 5
    u_opt5 = fmincon(@(uMV5) cost5(t(j), uMV5, u, p, s, x, 5, output, response), u_opt5,...
         [], [], [], [], p.MV_min_RP, p.MV_max_RP, [], options); 
    MV5(end+1) = u_opt5(1);
    output.MV5 = griddedInterpolant(t(1:j), MV5, 'previous');
    sol5 = odextend(sol5, @(t,x) FridgePlantsODEs5(s, p, x, u, t, output, 5,response), t(j+1));
    response.RP5(j,:) = deval(sol5, t(j));
    response.RP = [response.RP1 response.RP2 response.RP3 response.RP4 response.RP5];

    % Chill Dam
    u_opt_CD = fmincon(@(uMV) costCD(t(j), uMV, u, p, s, x, response, output), u_opt_CD, [], [], [], [],...
               p.MV_min_CD, p.MV_max_CD, [], options);
    MV_CD(end+1) = u_opt_CD(1); 
    output.MV_CD = griddedInterpolant(t(1:j), MV_CD, 'previous');
    solCD = odextend(solCD, @(t,x) ChillDamODEs(s, p, x, u, t, output, response), t(j+1)); 
    response.CD(:,j) = deval(solCD, t(j));

    T_UG(end+1) = v.T_outRPtot(end).*(v.F_outRPtot > 0.005) + x.T_CD(end).*(v.F_outRPtot <= 0.005);

        % if mod(j,(p.M_Econ*p.Ts)) == 0
        %     % Scheduling Optimisation (Economic Optimisation)
        %     GA_options = optimoptions('ga','Display','off', 'MaxGenerations', 1);
        %     u_opt_Econ = ga(@(uMV) costEcon(t(j), uMV, u, p, s, x, output), p.N_Econ,...
        %          [], [], [], [], p.MV_min_Econ, p.MV_max_Econ, ...
        %          @(uMV) mycon(t(j), uMV, x, output, s, u, z, p, response, t, ...
        %          solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, ...
        %          u_opt_SD, u_opt_PT, u_opt1, u_opt2, u_opt3, u_opt4, u_opt5, u_opt_CD),...
        %          1:p.N_Econ, GA_options);   
        %     output.MV_Econ(end+1) = u_opt_Econ(1);
        % end

        if mod(j,(p.N_Econ*p.Ts)) == 0
            RP_Cost = p.whole_week_tariff_vec_minutes(j+1);
            if RP_Cost < 0.9
                RP_Options = 5:-1:0;
                for plants = RP_Options
                    o = TestIfOver(t(j), x, output, s, u, z, p, response, t,...
                        solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, plants, u_opt_CD)
                    if o <= 100
                        output.MV_Econ(end+1) = plants;
                        break;
                    end
                end
            else
                RP_Options = 0:1:5;
                for plants = RP_Options
                    d = TestIfDry(t(j), x, output, s, u, z, p, response, t,...
                        solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, plants, u_opt_CD)
                    if d >= 20
                        output.MV_Econ(end+1) = plants;
                        break;
                    end
                end
            end
        end


        % if mod(j,(p.N_Econ*p.Ts)) == 0
        %     RP_Cost = p.whole_week_tariff_vec_minutes(j+1);
        %     if RP_Cost < 0.5
        %         RP_Options = 5:-1:0;
        %         for plants = RP_Options
        %             o = TestIfOver(t(j), x, output, s, u, z, p, response, t,...
        %                 solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, plants, u_opt_CD)
        %             if o <= 99
        %                 output.MV_Econ(end+1) = plants;
        %                 break;
        %             end
        %         end
        %     else if RP_Cost > 0.5 && RP_Cost < 0.9
        %             RP_Options = 1:1:3;
        %             for plants = RP_Options
        %                 d = TestIfDry(t(j), x, output, s, u, z, p, response, t,...
        %                     solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, plants, u_opt_CD)
        %                 if d >= 15
        %                     output.MV_Econ(end+1) = plants;
        %                     break;
        %                 end
        %             end
        %     else
        %         RP_Options = 0:1:5;
        %         for plants = RP_Options
        %             d = TestIfDry(t(j), x, output, s, u, z, p, response, t,...
        %                 solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, plants, u_opt_CD)
        %             if d >= 15
        %                 output.MV_Econ(end+1) = plants;
        %                 break;
        %             end
        %         end
        %     end
        %     end
        % end

    % saved.SP_L_SD(j,:) = output.MV_Econ(1,end); 
    % saved.SP_L_PT(j,:) = output.MV_Econ(2,end);
    % saved.SP_L_CD(j,:) = output.MV_Econ(3,end);
    % saved.SP_T_RP(j,:) = output.MV_Econ(4,end);   
    % saved.SP_T_CD(j,:) = output.MV_Econ(5,end);
    % saved.SP_no_RP(j,:) = output.MV_Econ(6,end);

    saved.SP_L_SD(j,:) = p.L_SD_SP_SS; 
    saved.SP_L_PT(j,:) = p.L_PT_SP_SS;
    saved.SP_L_CD(j,:) = p.L_CD_SP_SS;
    saved.SP_T_RP(j,:) = p.T_RP_SP_SS;   
    saved.SP_T_CD(j,:) = p.T_CD_SP_SS;
    saved.SP_no_RP(j,:) = output.MV_Econ(end);

    fprintf('%d\n',j) 

    end
%% Calculate Fridge Plant Intermediates
v.m_inRPtot = output.MV_PT(t(1:p.loop(end))) + output.MV_CD(t(1:p.loop(end)));
v.n         = saved.SP_no_RP; 
u_s = [1*(v.n>=1) 1*(v.n>1) 1*(v.n>2) 1*(v.n>3) 1*(v.n>4)];
v.n = sum(u_s,2);
v.F_outRP = (output.MV_CD(t(1:p.loop(end))) + output.MV_PT(t(1:p.loop(end))))./(v.n+0.0001).*u_s;
% if v.n > 0
%     v.F_outRP = (output.MV_CD(t(1:p.loop(end))) + output.MV_PT(t(1:p.loop(end))))./(v.n+0.001).*u_s; 
% else
%     v.F_outRP = u_s+0.001;
% end
v.F_Rec  = [output.MV1(t(1:p.loop(end))) output.MV2(t(1:p.loop(end))) output.MV3(t(1:p.loop(end))) output.MV4(t(1:p.loop(end))) output.MV5(t(1:p.loop(end)))];
v.F_inRP = v.F_outRP + v.F_Rec;       
v.T_inRPtot = ((output.MV_CD(t(1:p.loop(end))).*x.T_CD') + ...
              (output.MV_PT(t(1:p.loop(end))).*x.T_PT'))./...
              (output.MV_CD(t(1:p.loop(end))) + output.MV_PT(t(1:p.loop(end))));
v.T_RP   = [x.T1; x.T2; x.T3; x.T4; x.T5]'; 
v.T_inRP = ((v.F_Rec.*v.T_RP) + ...
           (v.F_outRP.*v.T_inRPtot))...
           ./ (v.F_Rec + v.F_outRP);                            
v.F_outRPtot = (v.F_outRP(:,1) + v.F_outRP(:,2) + v.F_outRP(:,3) + v.F_outRP(:,4) + v.F_outRP(:,5));
v.T_outRPtot = ((v.F_outRP(:,1).*x.T1') + (v.F_outRP(:,2).*x.T2') + (v.F_outRP(:,3).*x.T3') + (v.F_outRP(:,4).*x.T4')...
    + (v.F_outRP(:,5).*x.T5')) ./ v.F_outRPtot;
v.T_inRPtot = ((output.MV_CD(t(1:p.loop(end))).*response.CD(2,:)') + ...
              (output.MV_PT(t(1:p.loop(end))).*response.PT(2,:)'))./...
              (output.MV_CD(t(1:p.loop(end))) + output.MV_PT(t(1:p.loop(end))));
v.F_inCD = v.F_outRPtot - x.F_UG';
v.F_outCD =  output.MV_CD(t(1:p.loop(end))) + x.F_Ice';
v.F_inSysTot = output.MV_PT(t(1:p.loop(end)));
v.F_outSysTot = x.F_Ice' + x.F_UG';
v.T_UG = T_UG;

%% Plot MPC Results

% Square Dam
saved.SP_SD = saved.SP_L_SD;
lower_limit_SD = ones(1,length(p.loop)+1).*p.PV_min_SD(1);
upper_limit_SD = ones(1,length(p.loop)+1).*p.PV_max_SD(1);

figure(1)
ax1 = subplot(3,1,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, response.SD,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP_SD', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.L_SD, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.L_SD, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit_SD, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit_SD, '--k', 'MarkerSize', 3);
hold on
ylim([0 100]);
ylabel('Liquid level in SD (%)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax2 = subplot(3,1,2);
hold on
title('MV');
plot(t(1:p.loop(end))/p.Ts, output.MV_SD.Values(1:p.loop(end)),'b','LineWidth',1.5);
hold on
hold off
ax = gca;
ax.FontSize = p.font_size;
ylabel('F_o_u_t_S_D (L/S)'); xlabel('Time (min)');

ax3 = subplot(3,1,3);
hold on
title('DV');
plot(t(1:p.loop(end))/p.Ts, u.F_in_generatedSD(t(1:p.loop(end))),'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.F_in_SD, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.F_in_SD, '.k', 'MarkerSize', 5);
hold on
ylim([0 800]);
ylabel('F_i_n_S_D (L/S)'); xlabel('Time (min)');
legend('Ground Truth', 'State Estimate', 'Measurement');
ax = gca;
ax.FontSize = p.font_size;
hold off

linkaxes([ax1,ax2,ax3],'x');

% Pre-cooling Towers
saved.SP_PT = saved.SP_L_PT;
lower_limit_PT = ones(1,length(p.loop)+1).*p.PV_min_PT(1);
upper_limit_PT = ones(1,length(p.loop)+1).*p.PV_max_PT(1);

figure(2)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, response.PT(1,:),'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP_PT', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.L_PT, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.L_PT, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit_PT, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit_PT, '--k', 'MarkerSize', 3);
hold on
ylim([0 110]);
ylabel('Liquid level in PT (%)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:p.loop(end))/p.Ts, output.MV_PT.Values(1:p.loop(end)),'b','LineWidth',1.5);
hold on
ylim([0 900]);
ax = gca;
ax.FontSize = p.font_size;
hold off
ylabel('F_o_u_t_P_T (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
title('Non-temperature DVs');
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.F_AntiSurge_generated(t(1:p.loop(end))),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.F_AntiSurge, '--c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.F_AntiSurge, '.k', 'MarkerSize', 5);
hold on
ylabel('F_A_n_t_i_s_u_r_g_e (L/S)'); xlabel('Time (min)');
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.H_Env_generated(t(1:p.loop(end))),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.H_Env, '--','color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.H_Env, '.k', 'MarkerSize', 5);
hold on
ylabel('Humidity (%)'); xlabel('Time (min)');
legend('F_A_n_t_i_s_u_r_g_e Ground Truth', 'F_A_n_t_i_s_u_r_g_e State Estimate', 'F_A_n_t_i_s_u_r_g_e Measurement',...
       'H_E_n_v Ground Truth', 'H_E_n_v State Estimate', 'H_E_n_v Measurement');
ax = gca;
ax.FontSize = p.font_size;
hold off


ax4 = subplot(2,2,4);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840]);
title('Disturbance Temperatures and Outlet Temperature');
yyaxis left
hold on
plot(t(1:p.loop(end))/p.Ts, response.PT(2,:),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T_PT, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T_PT, '.k', 'MarkerSize', 5);
hold on
ylabel('T_P_T (^oC)'); xlabel('Time (min)');
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.T_AntiSurge_generated(t(1:p.loop(end))),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T_AntiSurge, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T_AntiSurge, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, u.T_out_generatedSD(t(1:p.loop(end))),'color',[0.9, 0.4, 0.09], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T_outSD, '--', 'color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T_outSD, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, u.T_amb_generated(t(1:p.loop(end))), 'color',[1, 0, 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T_amb,'--', 'color', [0.8, 0.05, 0.3], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T_amb, '.k', 'MarkerSize', 5);
hold on
ylabel('Disturbance Temperatures (^oC)'); xlabel('Time (min)');
legend('T_P_T Ground Truth','T_P_T State Estimate', 'T_P_T Measurement',...
       'T_A_n_t_i_s_u_r_g_e Ground Truth', 'T_A_n_t_i_s_u_r_g_e State Estimate', 'T_A_n_t_i_s_u_r_g_e Measurement',...
       'T_o_u_t_S_D Ground Truth', 'T_o_u_t_S_D State Estimate', 'T_o_u_t_S_D Measurement',...
       'T_a_m_b Ground Truth', 'T_a_m_b State Estimate', 'T_a_m_b Measurement');
ax = gca;
ax.FontSize = p.font_size;
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

% Fridge Plants
saved.SP_RP = saved.SP_T_RP;
lower_limit_RP = ones(1,length(p.loop)+1).*p.PV_min_RP(1);
upper_limit_RP = ones(1,length(p.loop)+1).*p.PV_max_RP(1);

F_Rec1 = output.MV1.Values';
T_in_ground1 = ((F_Rec1.*response.RP1) + ...
               (v.F_outRP(:,1).*v.T_inRPtot))...
              ./ (F_Rec1 + v.F_outRP(:,1)); 
F_Rec2 = output.MV2.Values';
T_in_ground2 = ((F_Rec2.*response.RP2) + ...
               (v.F_outRP(:,2).*v.T_inRPtot))...
              ./ (F_Rec2 + v.F_outRP(:,2)); 
F_Rec3 = output.MV3.Values';
T_in_ground3 = ((F_Rec3.*response.RP3) + ...
               (v.F_outRP(:,3).*v.T_inRPtot))...
              ./ (F_Rec3 + v.F_outRP(:,3)); 
F_Rec4 = output.MV4.Values';
T_in_ground4 = ((F_Rec4.*response.RP4) + ...
               (v.F_outRP(:,4).*v.T_inRPtot))...
              ./ (F_Rec4 + v.F_outRP(:,4)); 
F_Rec5 = output.MV5.Values';
T_in_ground5 = ((F_Rec5.*response.RP5) + ...
               (v.F_outRP(:,5).*v.T_inRPtot))...
              ./ (F_Rec5 + v.F_outRP(:,5)); 


% Fridge Plant 1
figure(3)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground1,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP_RP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin1, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin1, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit_RP, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit_RP, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_1 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:p.loop(end))/p.Ts, output.MV1.Values(1:p.loop(end)),'b','LineWidth',1.5);
hold on
ylim([0 150]);
ax = gca;
ax.FontSize = p.font_size;
hold off
ylabel('F_R_e_c_R_P_1 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:p.loop(end),1),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr1, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr1, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_1'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:p.loop(end),1),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb1, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb1, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_1'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_1 Ground Truth', 'Q_r_e_f_r_1 State Estimate', 'Q_r_e_f_r_1 Measurement', 'Q_a_m_b_1 Ground Truth', 'Q_a_m_b_1 State Estimate', 'Q_a_m_b_1 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperature');
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response.RP1,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T1, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T1, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_R_P_1 Ground Truth', 'T_R_P_1 State Estimate','T_R_P_1 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

% Fridge Plant 2
figure(4)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground2,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP_RP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin2, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin2, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit_RP, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit_RP, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_2 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:p.loop(end))/p.Ts, output.MV2.Values(1:p.loop(end)),'b','LineWidth',1.5);
hold on
ylim([0 150]);
ax = gca;
ax.FontSize = p.font_size;
hold off
ylabel('F_R_e_c_R_P_2 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:p.loop(end),2),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr2, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr2, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_2'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:p.loop(end),2),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb2, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb2, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_2'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_2 Ground Truth', 'Q_r_e_f_r_2 State Estimate', 'Q_r_e_f_r_2 Measurement', 'Q_a_m_b_2 Ground Truth', 'Q_a_m_b_2 State Estimate', 'Q_a_m_b_2 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperature');
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response.RP2,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T2, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T2, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_R_P_2 Ground Truth', 'T_R_P_2 State Estimate','T_R_P_2 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

% Fridge Plant 3
figure(5)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground3,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP_RP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin3, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin3, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit_RP, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit_RP, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_3 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:p.loop(end))/p.Ts, output.MV3.Values(1:p.loop(end)),'b','LineWidth',1.5);
hold on
ylim([0 150]);
ax = gca;
ax.FontSize = p.font_size;
hold off
ylabel('F_R_e_c_R_P_3 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:p.loop(end),3),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr3, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr3, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_3'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:p.loop(end),3),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb3, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb3, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_3'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_3 Ground Truth', 'Q_r_e_f_r_3 State Estimate', 'Q_r_e_f_r_3 Measurement', 'Q_a_m_b_3 Ground Truth', 'Q_a_m_b_3 State Estimate', 'Q_a_m_b_3 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperature');
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response.RP3,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T3, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T3, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_R_P_3 Ground Truth', 'T_R_P_3 State Estimate','T_R_P_3 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

% Fridge Plant 4
figure(6)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground4,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP_RP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin4, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin4, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit_RP, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit_RP, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_4 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:p.loop(end))/p.Ts, output.MV4.Values(1:p.loop(end)),'b','LineWidth',1.5);
hold on
ylim([0 150]);
ax = gca;
ax.FontSize = p.font_size;
hold off
ylabel('F_R_e_c_R_P_4 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:p.loop(end),4),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr4, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr4, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_4'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:p.loop(end),4),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb4, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb4, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_4'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_4 Ground Truth', 'Q_r_e_f_r_4 State Estimate', 'Q_r_e_f_r_4 Measurement', 'Q_a_m_b_4 Ground Truth', 'Q_a_m_b_4 State Estimate', 'Q_a_m_b_4 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperature');
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response.RP4,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T4, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T4, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_R_P_4 Ground Truth', 'T_R_P_4 State Estimate','T_R_P_4 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

% Fridge Plant 5
figure(7)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground5,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP_RP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin5, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin5, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit_RP, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit_RP, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_5 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:p.loop(end))/p.Ts, output.MV5.Values(1:p.loop(end)),'b','LineWidth',1.5);
hold on
ylim([0 150]);
ax = gca;
ax.FontSize = p.font_size;
hold off
ylabel('F_R_e_c_R_P_5 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:p.loop(end),5),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr5, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr5, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_5'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:p.loop(end),5),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb5, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb5, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_5'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_5 Ground Truth', 'Q_r_e_f_r_5 State Estimate', 'Q_r_e_f_r_5 Measurement', 'Q_a_m_b_5 Ground Truth', 'Q_a_m_b_5 State Estimate', 'Q_a_m_b_5 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperature');
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response.RP5,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T5, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T5, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_R_P_5 Ground Truth', 'T_R_P_5 State Estimate','T_R_P_5 Measurement', 'Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');


% Chill Dam
lower_limit_L_CD = ones(1,length(p.loop)+1).*p.PV_min_L_CD(1);
upper_limit_L_CD = ones(1,length(p.loop)+1).*p.PV_max_L_CD(1);

lower_limit_T_CD = ones(1,length(p.loop)+1).*p.PV_min_T_CD(1);
upper_limit_T_CD = ones(1,length(p.loop)+1).*p.PV_max_T_CD(1);


figure(8)
ax1 = subplot(2,2,1);
hold on
title('PVs');
colororder([0 0.5 0; 0 0.4470 0.7410])
yyaxis left
hold on
plot(t(1:p.loop(end))/p.Ts, response.CD(1,:),'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP_L_CD', 'r', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, x.L_CD, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.L_CD, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit_L_CD, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit_L_CD, '--k', 'MarkerSize', 3);
hold on
ylim([0 110]);
ylabel('Level_C_D (%)'); xlabel('Time (min)');
yyaxis right
hold on
plot(t(1:p.loop(end))/p.Ts, v.T_inRPtot,'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP_T_CD', 'r', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit_T_CD, '--b', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit_T_CD, '--b', 'MarkerSize', 3);
ylim([8 14]);
ylabel('T_i_n_R_P_t_o_t (^oC)'); xlabel('Time (min)');
hold on
legend('L_C_D Ground Truth','L_C_D SP', 'L_C_D State Estimate', 'L_C_D Measurement', 'L_C_D Lower Limit',...
       'L_C_D Upper Limit', 'T_i_n_R_P_t_o_t','T_i_n_R_P_t_o_t SP', 'T_i_n_R_P_t_o_t Lower Limit',...
       'T_i_n_R_P_t_o_t Upper Limit','Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:p.loop(end))/p.Ts, output.MV_CD.Values(1:p.loop(end)),'b','LineWidth',1.5);
hold on
ylim([0 750]);
ax = gca;
ax.FontSize = p.font_size;
hold off
ylabel('F_o_u_t_C_D (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
title('DVs');
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.F_Ice_generated(t(1:p.loop(end))),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.F_Ice, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.F_Ice, '.k', 'MarkerSize', 5);
hold on
ylim([0 1000]);
ylabel('F_I_c_e (L/S)'); xlabel('Time (min)');
legend('F_I_c_e Ground Truth', 'F_I_c_e State Estimate', 'F_I_c_e Measurement','Location', 'Best');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.F_UG_generated(t(1:p.loop(end))), 'color', [0.6350 0.0780 0.1840], 'LineWidth',1.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.F_UG, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.F_UG, '.k', 'MarkerSize', 5);
ylabel('F_U_G (L/S)'); xlabel('Time (min)');
legend('F_I_c_e Ground Truth','F_I_c_e State Estimate', 'F_I_c_e Measurement','F_U_G Ground Truth', 'F_U_G State Estimate', 'F_U_G Measurement','Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax4 = subplot(2,2,4);
hold on
title('Outlet Temperature');
plot(t(1:p.loop(end))/p.Ts, response.CD(2,:),'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T_CD, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T_CD, '.k', 'MarkerSize', 5);
hold on
ylim([0 15]);
ylabel('T_C_D (^oC)'); xlabel('Time (min)');
legend('Ground Truth', 'State Estimate', 'Measurement','Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

% Underground
figure(9)
ax1 = subplot(2,1,1);
hold on
title('Underground Flowrate');
plot(t(1:p.loop(end))/p.Ts, u.F_UG_generated(t(1:p.loop(end))),'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.F_UG, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.F_UG, '.k', 'MarkerSize', 5);
ylabel('F_U_G (L/s)'); xlabel('Time (min)');
legend('Ground Truth', 'State Estimate', 'Measurement','Location', 'Best');
ax = gca;
ax.FontSize = p.font_size;
hold off

ax2 = subplot(2,1,2);
hold on
title('Underground Temperature');
plot(t(1:p.loop(end))/p.Ts, v.T_UG,'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
ylabel('T_U_G (^oC)'); xlabel('Time (min)');
ax = gca;
ax.FontSize = p.font_size;
hold off
linkaxes([ax1,ax2],'x');

%% Cost Calculation and Plotting

Cost_RPs = p.whole_week_tariff_vec_minutes(1:p.Ts:length(p.loop)+1)'.*v.n(1:p.Ts:end).*p.pull_full_efficiency;
Total_Cost = sum(Cost_RPs);
Extrap_Cost_Month = Total_Cost*3.5*4;     % Cost for the month based on the simulation
Extrap_Cost_Year  = Extrap_Cost_Month*13; % Total cost for the year based on the simulation
table(Total_Cost, Extrap_Cost_Month, Extrap_Cost_Year)
format long;
figure(10)
colororder(['b'; 'k'])
yyaxis left
%plot(t(1:p.loop(end))/(p.Ts*60), Cost_RPs, 'LineWidth',1.5);
plot(t(1:p.loop(end))/(p.Ts*60), x.L_CD, 'LineWidth',1.5)
ax = gca;
ax.YAxis(1).Exponent = 0;
ax.FontSize = p.font_size;
hold on
ylabel('Chill Dam Level (%)'); xlabel('Time (h)');
hold on
yyaxis right
plot(t(1:p.loop(end))/(p.Ts*60), v.n, 'LineWidth',2.5);
ylim([0 5.5]);
ax = gca;
set(gca,'ytick', 0:5);
ax.FontSize = p.font_size;
ylabel('Number of RPs ON');
hold on
%tariff_change = findchangepts(p.whole_week_tariff_vec_minutes(1:length(p.loop)+2), MaxNumChanges=(length(p.loop)+2)/p.Ts/4)/p.Ts;
patch_y = [0 0 5.5 5.5];

patch_x_low1 = [0 6 6 0];
patch_x_low2 = [22 30 30 22];
patch_x_low3 = [46 48 48 46];
patch_x_med1 = [9 17 17 9];
patch_x_med2 = [19 22 22 19];
patch_x_med3 = [33 41 41 33];
patch_x_med4 = [43 46 46 43];
patch_x_hi1  = [6 9 9 6];
patch_x_hi2  = [17 19 19 17];
patch_x_hi3  = [30 33 33 30];
patch_x_hi4  = [41 43 43 41];

patch(patch_x_low1, patch_y, 'green','FaceAlpha', 0.1);
patch(patch_x_low2, patch_y, 'green','FaceAlpha', 0.1);
patch(patch_x_low3, patch_y, 'green','FaceAlpha', 0.1);

patch(patch_x_med1, patch_y, 'yellow','FaceAlpha', 0.1);
patch(patch_x_med2, patch_y, 'yellow','FaceAlpha', 0.1);
patch(patch_x_med3, patch_y, 'yellow','FaceAlpha', 0.1);
patch(patch_x_med4, patch_y, 'yellow','FaceAlpha', 0.1);

patch(patch_x_hi1, patch_y, 'red','FaceAlpha', 0.1);
patch(patch_x_hi2, patch_y, 'red','FaceAlpha', 0.1);
patch(patch_x_hi3, patch_y, 'red','FaceAlpha', 0.1);
patch(patch_x_hi4, patch_y, 'red','FaceAlpha', 0.1);

hold off

