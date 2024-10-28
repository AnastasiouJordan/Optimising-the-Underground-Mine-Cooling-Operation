%% System of ODEs: Refrigeration Plants
%  Jordan Anastasiou, 2023-11
%  This code is for the refrigeration plants where
%  each plant is modelled as a system, j.

clc
clear
clf

%% Define exogeneous inputs
load SavedInterpolantsRP.mat

%% Define Process Parameters
p.rho_Water = 1000;      % kg/m3,  Density of water
p.m_RPj     = 10000;     % kg,     Mass held by each fridge plant
                         %         Differs from v.m_RPj which is the mass
                         %         flowrate through each plant respectively
p.h_0       = 0.10186;   % kJ/kg,  Reference specific enthalpy
p.T_0       = 0.01;      % oC,     Reference temperature
p.C_p       = 4.1831;    % kJ/kgC, Heat capacity of water

% % Initial guesses for unknown parameters UA_RPj and UA_amb
p.UA_RP1  = 13000;       % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 1
p.UA_RP2  = 10000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 2
p.UA_RP3  = 15000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 3
p.UA_RP4  = 32500;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 4
p.UA_RP5  = 23000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 5
p.UA_RP   = [p.UA_RP1, p.UA_RP2, p.UA_RP3, p.UA_RP4, p.UA_RP5];

p.UA_amb1 = -100;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 1
p.UA_amb2 = -500;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 2
p.UA_amb3 = 200;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 3
p.UA_amb4 = -300;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 4
p.UA_amb5 = -500;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 5
p.UA_amb  = [p.UA_amb1, p.UA_amb2, p.UA_amb3, p.UA_amb4, p.UA_amb5];


% Initial guesses for unknown parameters UA_RPj and UA_amb
% p.UA_RP1  = 10000;       % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 1
% p.UA_RP2  = 7500;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 2
% p.UA_RP3  = 7500;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 3
% p.UA_RP4  = 18000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 4
% p.UA_RP5  = 23000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 5
% p.UA_RP   = [p.UA_RP1, p.UA_RP2, p.UA_RP3, p.UA_RP4, p.UA_RP5];
% 
% p.UA_amb1 = 100;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 1
% p.UA_amb2 = 275;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 2
% p.UA_amb3 = 275;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 3
% p.UA_amb4 = 100;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 4
% p.UA_amb5 = 100;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 5
% p.UA_amb  = [p.UA_amb1, p.UA_amb2, p.UA_amb3, p.UA_amb4, p.UA_amb5];

p.regressedparameterfields = {'UA_RP', 'UA_amb'}; % Parameter structure field names

pUAvec = S2V(p, p.regressedparameterfields); % Convert the unknown parameters to a vector of unknown parameters

%% Load AR Model Data
% Load the AR model data, providing stochastic sets of data
% representing the trends in the measurement data, with
% calculated variances and constants.
[p,u] = ARModels(t,u,p);

%% Define Kalman Filter parameters
p.xhat_0 = [4; 1; -0.05; 12];         % Intial state estimates: T_RPj, T_refrj, T_amb, T_inRPj
p.H      = eye(length(p.xhat_0));  % Observation matrix 
p.c_K = [0 0 0 0 0; p.c_Q_refr; p.c_Q_amb; p.c_T_in]; % Constants from AR model
p.w_K = [p.w_T ; p.w_Q_refr; p.w_Q_amb; p.w_T_in];    % Variance from AR model
p.T_sensor_noise   = 0.01;         % Std dev of measurement noise, level sensor error
p.Q_lumped_noise   = 0.01;         % Std dev of measurement noise, inlet stream flowmeter error
p.R = diag([p.T_sensor_noise^2 ;...
            p.Q_lumped_noise^2;...
            p.Q_lumped_noise^2 ;...
            p.T_sensor_noise^2 ]);  % Measurement noise covariance matrix (measurement error squared) 
p.Q1 = diag([p.w_T(1); 0.002;  1E-12; 1E-12]); % Process noise covariance matrix based on variance from AR model for RP1
p.Q2 = diag([p.w_T(2); 0.02;  1E-12; 1E-12]); % Process noise covariance matrix based on variance from AR model for RP2
p.Q3 = diag([p.w_T(3); 0.02;  1E-12; 1E-12]); % Process noise covariance matrix based on variance from AR model for RP3
p.Q4 = diag([p.w_T(4); 0.02;  1E-12; 1E-12]); % Process noise covariance matrix based on variance from AR model for RP4
p.Q5 = diag([p.w_T(5); 0.02;  1E-12; 1E-12]); % Process noise covariance matrix based on variance from AR model for RP5

%% Define MPC parameters
p.T_SS       = 13;                          % oC, Steady state inlet temperature of the fridge plants (initial PV value)
p.N          = 3;                           % ~, Number of samples of the control input/prediction nodes
p.Ts         = 60;                          % s, Sampling period, the frequency at which a new control input is determined
p.Stp        = length(t);                   % ~, Number of steps in simulation
p.TL         = (p.Ts*p.Stp) - p.Ts;         % s, Total time or Time Limit (sampling period x total time)
p.loop       = 2:1:400;                      % ~, Simulation points for the loop (mostly for testing)
p.F_RecSS    = 50;                          % L/s, Steady state recycle flowrate (initial MV value)
p.uvec_init  = p.F_RecSS*ones(1,p.N);       % L/s, Initial points (sequence guess)
p.SP         = 11*ones(1,p.N);              % oC, Initial SP for the outlet temperature of the fridge plants (initial SP value)
p.SP_changes = 2296;                        % ~, Number of SP changes (2296 - every 22 min. 1335 - every 38 min)
p.SP_min     = 10;                          % %, Lowest SP for the outlet temperature of the fridge plants
p.SP_max     = 12;                          % %, Highest SP for the outlet temperature of the fridge plants
p.SP_samples = p.SP_min...
               + (p.SP_max - p.SP_min)...
               * rand(p.SP_changes, 1);     % Sample SP changes
p.SP_times   = (0:p.TL/p.SP_changes:p.TL)'; % Times at which the SP should change
p.MV_min     = 0*ones(1,p.N);               % L/s, Minimum MV limit (ensure that there is always enough flow for RPs)
p.MV_max     = 100*ones(1,p.N);             % L/s, Maximum MV limit
p.PV_min     = 10*ones(1,p.N);              % %, Minimum PV limit (minimum temperature)
p.PV_max     = 12*ones(1,p.N);              % %, Maximum PV limit (maximum temperature)
p.Q_Weight   = 5;                           % SP weight
p.R_Weight   = 0.1;                         % MV weight

%% Define state structure and initial conditions

s.statefields    = {'T1',      'T2',      'T3',      'T4',      'T5'};
s.MPCstatefields = {'T1',      'T2',      'T3',      'T4',      'T5';...
                    'Qrefr1',  'Qrefr2',  'Qrefr3',  'Qrefr4',  'Qrefr5';...
                    'Qamb1',   'Qamb2',   'Qamb3',   'Qamb4',   'Qamb5';...
                    'Tin1',    'Tin2',    'Tin3',    'Tin4',    'Tin5'}; % Field names for each state in MPC
s.KFstatefields  = {'T1',      'T2',      'T3',      'T4',      'T5';...
                    'Qrefr1',  'Qrefr2',  'Qrefr3',  'Qrefr4',  'Qrefr5';...
                    'Qamb1',   'Qamb2',   'Qamb3',   'Qamb4',   'Qamb5';...
                    'Tin1',    'Tin2',    'Tin3',    'Tin4',    'Tin5';...
                    'P1',      'P2',      'P3',      'P4',      'P5'};   % Field names for each state in KF

% Initial values for the fridge plant outlet temperatures
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

x0.Tin1 = u.T_in_generated.Values(1,1);  
x0.Tin2 = u.T_in_generated.Values(1,2);
x0.Tin3 = u.T_in_generated.Values(1,3);
x0.Tin4 = u.T_in_generated.Values(1,4);
x0.Tin5 = u.T_in_generated.Values(1,5);

x0.P    = 0.1;         

save FridgePlants.mat s x0 p u t

%%
clc
clear
load FridgePlants.mat

%% MPC Initialisation
% Function that generates the optimal sequence of
% actions given the currrent starting state and SP
% (initialise and run for first time step)

options = optimoptions('fmincon','Display','off','EnableFeasibilityMode',true);

% Initialise the solution structure for each fridge plant
sol1.y = [x0.T1];
sol2.y = [x0.T2];
sol3.y = [x0.T3];
sol4.y = [x0.T4];
sol5.y = [x0.T5];

% Initialise measurement structure
z.T1     = [];     z.T2 = [];     z.T3 = [];     z.T4 = [];     z.T5 = [];
z.Qrefr1 = []; z.Qrefr2 = []; z.Qrefr3 = []; z.Qrefr4 = []; z.Qrefr5 = []; 
z.Qamb1  = [];  z.Qamb2 = [];  z.Qamb3 = [];  z.Qamb4 = [];  z.Qamb5 = []; 
z.Tin1   = [];   z.Tin2 = [];   z.Tin3 = [];   z.Tin4 = [];   z.Tin5 = [];


v.T_inRP = [p.T_SS p.T_SS p.T_SS p.T_SS p.T_SS];
z    = Meas(sol1, sol2, sol3, sol4, sol5, u, 0, p, z, v); 
x    = z;   % Set the state estimates to initially be equal to the measurements

% Set the initial covariance to be used in the Kalman Filter
x.P1 = 0.1; x.P2 = 0.1; x.P3 = 0.1; x.P4 = 0.1; x.P5 = 0.1;

% Fridge Plant 1
% Calculate the optimal MV movement
u_opt1 = fmincon(@(uMV) cost1(t(1), uMV, u, p, s, x, 1), p.uvec_init,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
MV1    = u_opt1(1);
output.MV1 = @(t) MV1; 
% Simulate the plant's response to the provided MV
sol1 = ode45(@(t,x) FridgePlantsODEs1(s,p,x,u,t,output,1), [t(1) t(2)], sol1.y); % Calculate the ground truth for the first time interval
response1(1,:) = deval(sol1, t(1));

% Fridge Plant 2
% Calculate the optimal MV movement
u_opt2 = fmincon(@(uMV) cost2(t(1), uMV, u, p, s, x, 2), p.uvec_init,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
MV2    = u_opt2(1);
output.MV2 = @(t) MV2; 
% Simulate the plant's response to the provided MV
sol2 = ode45(@(t,x) FridgePlantsODEs2(s,p,x,u,t,output,2), [t(1) t(2)], sol2.y); % Calculate the ground truth for the first time interval
response2(1,:) = deval(sol2, t(1));

% Fridge Plant 3
% Calculate the optimal MV movement
u_opt3 = fmincon(@(uMV) cost3(t(1), uMV, u, p, s, x, 3), p.uvec_init,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
MV3    = u_opt3(1);
output.MV3 = @(t) MV3; 
% Simulate the plant's response to the provided MV
sol3 = ode45(@(t,x) FridgePlantsODEs3(s,p,x,u,t,output,3), [t(1) t(2)], sol3.y); % Calculate the ground truth for the first time interval
response3(1,:) = deval(sol3, t(1));

% Fridge Plant 4
% Calculate the optimal MV movement
u_opt4 = fmincon(@(uMV) cost4(t(1), uMV, u, p, s, x, 4), p.uvec_init,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
MV4    = u_opt4(1);
output.MV4 = @(t) MV4; 
% Simulate the plant's response to the provided MV
sol4 = ode45(@(t,x) FridgePlantsODEs4(s,p,x,u,t,output,4), [t(1) t(2)], sol4.y); % Calculate the ground truth for the first time interval
response4(1,:) = deval(sol4, t(1));

% Fridge Plant 5
% Calculate the optimal MV movement
u_opt5 = fmincon(@(uMV) cost5(t(1), uMV, u, p, s, x, 5), p.uvec_init,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
MV5    = u_opt5(1);
output.MV5 = @(t) MV5; 
% Simulate the plant's response to the provided MV
sol5 = ode45(@(t,x) FridgePlantsODEs5(s,p,x,u,t,output,5), [t(1) t(2)], sol5.y); % Calculate the ground truth for the first time interval
response5(1,:) = deval(sol5, t(1));

response = [];
response = [response1 response2 response3 response4 response5];
for b = 1:2
    saved.SP_init(b,:) = p.SP; % Save SPs
end

%% MPC Loop

for j = p.loop % Using shorter loop for testing to reduce run-time
        
    % Change SP and save it
    for sp = 1:1:size(p.SP_times,1)
        if p.SP_times(sp) == j*p.Ts
            p.SP = p.SP_samples(sp)*ones(1,p.N);
        else
            p.SP = p.SP;
        end
    end
    saved.SP_loop(j,:) = p.SP; 
    
    v = RPIntermediatesKF(x, u, p, t(j), output, s, response);
    z = Meas(sol1, sol2, sol3, sol4, sol5, u, t(j), p, z, v);
    x = KalmanFilterRP(x, s, output, u, z, [t(j-1) t(j)], p, response);

    % Fridge Plant 1
    % Calculate the optimal MV movement
    u_opt1 = fmincon(@(uMV1) cost1(t(j), uMV1, u, p, s, x, 1), u_opt1,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
    MV1(end+1) = u_opt1(1);
    output.MV1 = griddedInterpolant(t(1:j), MV1, 'previous');
    % Simulate the plant's response to the provided MV
    sol1 = odextend(sol1, @(t,x) FridgePlantsODEs1(s, p, x, u, t, output, 1), t(j+1)); % Ground truth (how the system is actually responding)
    response1(j,:) = deval(sol1, t(j));
   
    % Fridge Plant 2
    % Calculate the optimal MV movement
    u_opt2 = fmincon(@(uMV2) cost2(t(j), uMV2, u, p, s, x, 2), u_opt2,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
    MV2(end+1) = u_opt2(1);
    output.MV2 = griddedInterpolant(t(1:j), MV2, 'previous');
    % Simulate the plant's response to the provided MV
    sol2 = odextend(sol2, @(t,x) FridgePlantsODEs2(s, p, x, u, t, output, 2), t(j+1)); % Ground truth (how the system is actually responding)
    response2(j,:) = deval(sol2, t(j));

    % Fridge Plant 3
    % Calculate the optimal MV movement
    u_opt3 = fmincon(@(uMV3) cost3(t(j), uMV3, u, p, s, x, 3), u_opt3,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
    MV3(end+1) = u_opt3(1);
    output.MV3 = griddedInterpolant(t(1:j), MV3, 'previous');
    % Simulate the plant's response to the provided MV
    sol3 = odextend(sol3, @(t,x) FridgePlantsODEs3(s, p, x, u, t, output, 3), t(j+1)); % Ground truth (how the system is actually responding)
    response3(j,:) = deval(sol3, t(j));

    % Fridge Plant 4
    % Calculate the optimal MV movement
    u_opt4 = fmincon(@(uMV4) cost4(t(j), uMV4, u, p, s, x, 4), u_opt4,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
    MV4(end+1) = u_opt4(1);
    output.MV4 = griddedInterpolant(t(1:j), MV4, 'previous');
    % Simulate the plant's response to the provided MV
    sol4 = odextend(sol4, @(t,x) FridgePlantsODEs4(s, p, x, u, t, output, 4), t(j+1)); % Ground truth (how the system is actually responding)
    response4(j,:) = deval(sol4, t(j));

    % Fridge Plant 5
    % Calculate the optimal MV movement
    u_opt5 = fmincon(@(uMV5) cost5(t(j), uMV5, u, p, s, x, 5), u_opt5,...
         [], [], [], [],p.MV_min,p.MV_max, [], options); 
    MV5(end+1) = u_opt5(1);
    output.MV5 = griddedInterpolant(t(1:j), MV5, 'previous');
    % Simulate the plant's response to the provided MV
    sol5 = odextend(sol5, @(t,x) FridgePlantsODEs5(s, p, x, u, t, output, 5), t(j+1)); % Ground truth (how the system is actually responding)
    response5(j,:) = deval(sol5, t(j));
      
    response = [response1 response2 response3 response4 response5];
    fprintf('%d\n',j) 
end
v = RPIntermediates(x,u,p,t(1:p.loop(end)),output,s);

% Calculate the ground truth inlet temperature
F_Rec1 = output.MV1.Values';
T_in_ground1 = ((F_Rec1.*response1) + ...
               (v.F_outRP(:,1).*v.T_inRPtot))...
              ./ (F_Rec1 + v.F_outRP(:,1)); 
F_Rec2 = output.MV2.Values';
T_in_ground2 = ((F_Rec2.*response2) + ...
               (v.F_outRP(:,2).*v.T_inRPtot))...
              ./ (F_Rec2 + v.F_outRP(:,2)); 
F_Rec3 = output.MV3.Values';
T_in_ground3 = ((F_Rec3.*response3) + ...
               (v.F_outRP(:,3).*v.T_inRPtot))...
              ./ (F_Rec3 + v.F_outRP(:,3)); 
F_Rec4 = output.MV4.Values';
T_in_ground4 = ((F_Rec4.*response4) + ...
               (v.F_outRP(:,4).*v.T_inRPtot))...
              ./ (F_Rec4 + v.F_outRP(:,4)); 
F_Rec5 = output.MV5.Values';
T_in_ground5 = ((F_Rec5.*response5) + ...
               (v.F_outRP(:,5).*v.T_inRPtot))...
              ./ (F_Rec5 + v.F_outRP(:,5)); 

saved.SP(1:size(saved.SP_init,1)) = saved.SP_init(:,1);
saved.SP(size(saved.SP_init,1)+1:length(saved.SP_loop)+1) = saved.SP_loop(2:end,1);
saved.SP = saved.SP(1:end-1);

% Level Limits
lower_limit = ones(1,length(p.loop)+1).*p.PV_min(1);
upper_limit = ones(1,length(p.loop)+1).*p.PV_max(1);

%% Fridge Plant 1
figure(1)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground1,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin1, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin1, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_1 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:j)/p.Ts, output.MV1.Values(1:j),'b','LineWidth',1.5);
hold on
ylim([0 100]);
hold off
ylabel('F_R_e_c_R_P_1 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:j,1),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr1, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr1, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_1'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:j,1),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb1, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb1, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_1'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_1 Ground Truth', 'Q_r_e_f_r_1 State Estimate', 'Q_r_e_f_r_1 Measurement', 'Q_a_m_b_1 Ground Truth', 'Q_a_m_b_1 State Estimate', 'Q_a_m_b_1 Measurement', 'Location', 'Best');
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperatures');
plot(t(1:p.loop(end))/p.Ts, u.T_outCD.Values(1:j),'.','Color',[0.8500 0.3250 0.0980]); 
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response1,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T1, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T1, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_o_u_t_C_D', 'T_R_P_1 Ground Truth', 'T_R_P_1 State Estimate','T_R_P_1 Measurement', 'Location', 'Best');
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

%% Fridge Plant 2
figure(2)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground2,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin2, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin2, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_2 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:j)/p.Ts, output.MV2.Values(1:j),'b','LineWidth',1.5);
hold on
ylim([0 100]);
hold off
ylabel('F_R_e_c_R_P_2 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:j,2),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr2, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr2, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_2'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:j,2),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb2, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb2, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_2'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_2 Ground Truth', 'Q_r_e_f_r_2 State Estimate', 'Q_r_e_f_r_2 Measurement', 'Q_a_m_b_2 Ground Truth', 'Q_a_m_b_2 State Estimate', 'Q_a_m_b_2 Measurement', 'Location', 'Best');
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperatures');
plot(t(1:p.loop(end))/p.Ts, u.T_outCD.Values(1:j),'.','Color',[0.8500 0.3250 0.0980]); 
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response2,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T2, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T2, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_o_u_t_C_D', 'T_R_P_2 Ground Truth', 'T_R_P_2 State Estimate','T_R_P_2 Measurement', 'Location', 'Best');
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

%% Fridge Plant 3
figure(3)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground3,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin3, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin3, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_3 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:j)/p.Ts, output.MV3.Values(1:j),'b','LineWidth',1.5);
hold on
ylim([0 100]);
hold off
ylabel('F_R_e_c_R_P_3 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:j,3),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr3, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr3, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_3'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:j,3),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb3, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb3, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_3'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_3 Ground Truth', 'Q_r_e_f_r_3 State Estimate', 'Q_r_e_f_r_3 Measurement', 'Q_a_m_b_3 Ground Truth', 'Q_a_m_b_3 State Estimate', 'Q_a_m_b_3 Measurement', 'Location', 'Best');
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperatures');
plot(t(1:p.loop(end))/p.Ts, u.T_outCD.Values(1:j),'.','Color',[0.8500 0.3250 0.0980]); 
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response3,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T3, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T3, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_o_u_t_C_D', 'T_R_P_3 Ground Truth', 'T_R_P_3 State Estimate','T_R_P_3 Measurement', 'Location', 'Best');
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

%% Fridge Plant 4
figure(4)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground4,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin4, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin4, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_4 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:j)/p.Ts, output.MV4.Values(1:j),'b','LineWidth',1.5);
hold on
ylim([0 100]);
hold off
ylabel('F_R_e_c_R_P_4 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:j,4),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr4, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr4, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_4'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:j,4),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb4, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb4, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_4'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_4 Ground Truth', 'Q_r_e_f_r_4 State Estimate', 'Q_r_e_f_r_4 Measurement', 'Q_a_m_b_4 Ground Truth', 'Q_a_m_b_4 State Estimate', 'Q_a_m_b_4 Measurement', 'Location', 'Best');
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperatures');
plot(t(1:p.loop(end))/p.Ts, u.T_outCD.Values(1:j),'.','Color',[0.8500 0.3250 0.0980]); 
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response4,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T4, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T4, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_o_u_t_C_D', 'T_R_P_4 Ground Truth', 'T_R_P_4 State Estimate','T_R_P_4 Measurement', 'Location', 'Best');
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');

%% Fridge Plant 5
figure(5)
ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, T_in_ground5,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.Tin5, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Tin5, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit, '--k', 'MarkerSize', 3);
hold on
ylim([9 16]);
ylabel('T_i_n_R_P_5 (^oC)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:j)/p.Ts, output.MV5.Values(1:j),'b','LineWidth',1.5);
hold on
ylim([0 100]);
hold off
ylabel('F_R_e_c_R_P_5 (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
colororder([0 0.4470 0.7410; 0.6350 0.0780 0.1840])
title('DVs');
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.Q_refr_generated.Values(1:j,5),'color',[0 0.4470 0.7410], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qrefr5, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qrefr5, '.k', 'MarkerSize', 5);
hold on
ylim([0 4]);
ylabel('Q_r_e_f_r_5'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.Q_amb_generated.Values(1:j,5),'color',[0.6350 0.0780 0.1840], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.Qamb5, 'color', [0.8500 0.3250 0.0980], 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.Qamb5, '.k', 'MarkerSize', 5);
ylabel('Q_a_m_b_5'); xlabel('Time (min)');
hold on
legend('Q_r_e_f_r_5 Ground Truth', 'Q_r_e_f_r_5 State Estimate', 'Q_r_e_f_r_5 Measurement', 'Q_a_m_b_5 Ground Truth', 'Q_a_m_b_5 State Estimate', 'Q_a_m_b_5 Measurement', 'Location', 'Best');
hold off


ax4 = subplot(2,2,4);
hold on
title('Temperatures');
plot(t(1:p.loop(end))/p.Ts, u.T_outCD.Values(1:j),'.','Color',[0.8500 0.3250 0.0980]); 
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, response5,'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.T5, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.T5, '.k', 'MarkerSize', 5);
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_o_u_t_C_D', 'T_R_P_5 Ground Truth', 'T_R_P_5 State Estimate','T_R_P_5 Measurement', 'Location', 'Best');
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');
