%% System of ODEs: Pre-cooling Towers
%  Jordan Anastasiou, 2022-06
%  This code is for the pre-cooling towers with both towers
%  being grouped together to form one system over which the
%  mass and energy balances were performed. There is a flowrate
%  of water as well as a flowrate of air through the tower.
clc
clear
clf

%% Define exogeneous inputs
load SavedInterpolantsPT.mat

%% Define Process Parameters
p.regressedparameterfields = {'m_Air'}; 

p.rho_Water = 1000;      % kg/m3,  Density of water
p.rho_Air   = 1.225;     % kg/m3,  Density of air
p.H_outAir  = 100;       % %,      Humidity of air leaving towers
p.h_0       = 0.10186;   % kJ/kg,  Reference specific enthalpy
p.T_0       = 0.01;      % oC,     Reference temperature
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
p.C_p       = 4.1831;    % kJ/kgC, Heat capacity of water

p.height_PT = 3;         % m, Height of the Pre-cooling Tower basins
p.length_PT = 14.8;      % m, Length of the Pre-cooling Tower basins
p.width_PT  = 10.45*2;   % m, Width of the Pre-cooling Tower basins
p.area_PT   = p.length_PT*p.width_PT; % m2, Area of the Pre-cooling Tower basins
p.volume_PT = p.area_PT*p.height_PT;  % m3, Volume of the Pre-cooling Tower basins

p.m_PTAmax  = p.volume_PT/2*1000;     % kg, Maximum mass capacity of PT A basin
p.m_PTBmax  = p.volume_PT/2*1000;     % kg, Maximum mass capacity of PT B basin
p.m_PTmax   = p.m_PTAmax*2;

% Initial guess for m_Air
p.m_Air     = 800;       % kg/s, Estimated mass flowrate of air.
                         % This is based on a L/G ratio of 0.5 from 
                         % literature with the water flowrate through
                         % the towers having an average range of 
                         % 300 - 500 kg/s.
pmAirVec = S2V(p, p.regressedparameterfields); % Convert the unknown parameter to a vector

%% Load AR Model Data
% Load the AR model data, providing stochastic sets of data
% representing the trends in the measurement data, with
% calculated variances and constants.
[p,u] = ARModels(t,u,p);

%% Define Kalman Filter parameters
p.xhat_0 = [2.9; 400];              % Intial state estimates (Level in m, inlet flowrate in L/s)
p.H      = eye(length(p.xhat_0));   % Observation matrix 
% p.H = [1 0 0;...                  % Observation matrix with constraint (for constrained KF only)
%      0 1 0;...
%      0 0 1;...
%      1 -deltaT/2 deltaT/2];
p.c_K = [0 ; p.c_F_in];            % Constants from AR model
p.w_K = [p.w_L ; p.w_F_in];       % Variance from AR model
p.L_sensor_noise   = 0.01;          % Std dev of measurement noise, level sensor error
p.F_in_meter_noise = 0.01;          % Std dev of measurement noise, inlet stream flowmeter error
p.R = diag([p.L_sensor_noise^2 ;...
            p.F_in_meter_noise^2]); % Measurement noise covariance matrix (measurement error squared)
% p.Q = diag([p.w_K(1) ; p.w_K(2)]);% Process noise covariance matrix based on variance from AR model   
p.Q = diag([p.w_K(1); p.w_K(2)]);  % Process noise covariance matrix based on variance from AR model 
% p.A = [1 p.C_L;    
%        0 a.F_in];                 % Transition matrix
p.A = [0 p.C_L;    
       0 0];                        % Transition matrix

%% Define MPC Parameters
p.L_SS       = 80;                          % %, Steady state level in the square dam (initial PV value)
p.N          = 3;                           % ~, Number of samples of the control input/prediction nodes
p.Ts         = 60;                          % s, Sampling period, the frequency at which a new control input is determined
p.Stp        = length(t);                   % ~, Number of steps in simulation
p.TL         = (p.Ts*p.Stp) - p.Ts;         % s, Total time or Time Limit (sampling period x total time)
p.loop       = 2:1:500;                      % ~, Simulation points for the loop (mostly for testing)
p.F_outSS    = u.F_outSD(1);                % L/s, Steady state outlet flowrate (initial MV value)
p.uvec_init  = p.F_outSS*ones(1,p.N);       % L/s, Initial points (sequence guess)
p.SP         = 85*ones(1,p.N);              % %, Initial SP for the level in the square dam (initial SP value)
p.SP_changes = 2296;                        % ~, Number of SP changes (2296 - every 22 min. 1335 - every 38 min)
p.SP_min     = 75;                          % %, Lowest SP for the level in the Dam
p.SP_max     = 90;                          % %, Highest SP for the level in the Dam
p.SP_samples = p.SP_min...
               + (p.SP_max - p.SP_min)...
               * rand(p.SP_changes, 1);     % Sample SP changes
p.SP_times   = (0:p.TL/p.SP_changes:p.TL)'; % Times at which the SP should change
p.MV_min     = 300*ones(1,p.N);             % L/s, Minimum MV limit (ensure that there is always enough flow for RPs)
p.MV_max     = 800*ones(1,p.N);             % L/s, Maximum MV limit
p.PV_min     = 40*ones(1,p.N);              % %, Minimum PV limit (minimum level)
p.PV_max     = 95*ones(1,p.N);              % %, Maximum PV limit (maximum level)
p.Q_Weight   = 1;                           % SP weight
p.R_Weight   = 0.1;                         % MV weight
%% Define state structure and initial conditions

s.statefields    = {'L','h_PT'};              % Field names for each state  
s.MPCstatefields = {'L','F_in', 'h_PT'};      % Field names for each state in MPC
s.KFstatefields  = {'L','F_in', 'h_PT', 'P'}; % Field names for each state in KF

x0.L    = (u.L_PTA(0) + u.L_PTB(0))/2;     % %, Initial value level
x0.F_in = u.F_outSD(0) + u.F_AntiSurge(0); % L/s, Initial value for inlet flowrate

x0.h_PT = 80;  % kJ/kgK, Initial value for enthalpy of water leaving PTs
x0.P    = 0.1; % Initial state estimate covariance matrix


x0_vec  = S2V(x0, s.statefields);

%% Generate variable inlet flowrate
% In the ARIMA Models file, autoregressive data is produced for use in the
% KF. An additional set is created using the same model to be used in the 
% MPC. This has been named u.F_in_generated.

u.F_in_generated(t);

%% Simulate system of ODEs
% [~, x_vec] = ode45(@(t, x) PreCoolingTowersODEs(s, p, x, u, t), t, x0_vec);
% x = V2S(x_vec', s.statefields);
% v = PTIntermediates(x, u, p, t);
% save PreCoolingTowers.mat v u

%% Plot
% State variables. Plot dm_PT/dt and dh_PT/dt calculated using ODE.
font_size = 25;

% figure (1)
% title('Prediction Results with Assumed Parameter Values');
% subplot(2,1,1)
% plot(t/86400, u.T_PT(t), t/86400, v.T_PT);  
% legend('measured', 'predicted', 'FontSize', font_size);
% %text(3100000, 2, sprintf('mAir = %.2f kJ/Ks', p.m_Air));
% xlabel('Time (days)');
% ylabel('T_P_T (^oC)');
% ax = gca;
% ax.FontSize = font_size;
% subplot(2,1,2)
% plot(t/86400, u.L_PTA(t), t/86400, v.L_PT);
% legend('measured', 'predicted', 'FontSize', font_size);
% %text(3100000, 2, sprintf('mAir = %.2f kJ/Ks', p.m_Air));
% xlabel('Time (days)');
% ylabel('L_P_T (%)');
% ax = gca;
% ax.FontSize = font_size;


%% Regression

% options = optimoptions('lsqnonlin', 'StepTolerance', 1e-9,...
%                        'Algorithm','trust-region-reflective',...
%                        'FiniteDifferenceType','central',...
%                        'TypicalX', 650,...
%                        'Display', 'iter-detailed');
% 
% p_est    = lsqnonlin(@(pmAirVec) PTCalcError(pmAirVec, u, p, s, t), pmAirVec, 100, 900, options)
% 
% [E, x, v] = PTCalcError(p_est, u, p, s, t);
% 
% %Plot results using parameter estimate/regressed parameter
% figure (2)
% title('Prediction Results with Regressed Parameter Values');
% subplot(2,1,1)
% plot(t, u.T_PT(t), t, v.T_PT);  
% legend('measured', 'predicted');
% text(3100000, 2, sprintf('mAir = %.2f kJ/Ks', p.m_Air));
% xlabel('Time (s)');
% ylabel('T_P_T (^oC)');
% subplot(2,1,2)
% plot(t, u.L_PTA(t), t, v.L_PT);
% legend('measured', 'predicted');
% text(3100000, 2, sprintf('mAir = %.2f kJ/Ks', p.m_Air));
% xlabel('Time (s)');
% ylabel('L_P_T (%)');
% ax = gca;
% ax.FontSize = font_size;

%% Likelihood Profiles

% SSR_Level = []; % Sum of squared residuals for level
% SSR_Temp  = []; % Sum of squared residuals for temperature
% soln_space = 0.01:200:3000;
% for pmAirvec = soln_space
%     [E, x, v] = PTCalcError(pmAirvec, u, p, s, t);
%     E_2 = E.^2;
%     sum_E_2_Level = sum(E_2(:,1));
%     sum_E_2_Temp  = sum(E_2(:,2));
%     SSR_Level = [SSR_Level sum_E_2_Level];
%     SSR_Temp  = [SSR_Temp sum_E_2_Temp];
% end
% 
% figure(3)
% title('Likelihood Profile for parameter mAir and Level');
% LLRatio_Level = 2*log(SSR_Level/min(SSR_Level));
% plot(soln_space, LLRatio_Level);
% hold on
% yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
% hold off
% xlabel('m_A_i_r (m/s)');
% ylabel('Negative Log Likelihood Ratio');
% xlim([0 soln_space(end)]);
% x1 = interp1(LLRatio_Level, soln_space, 0);
% zero_point = find(soln_space == x1);
% x2 = interp1(LLRatio_Level(1:zero_point), soln_space(1:zero_point), 2.71);
% x3 = interp1(LLRatio_Level(zero_point:end), soln_space(zero_point:end), 2.71);
% xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
% xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
% xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
% ax = gca;
% ax.FontSize = font_size;
% 
% figure(4)
% LLRatio_Temp = 2*log(SSR_Temp/min(SSR_Temp));
% plot(soln_space, LLRatio_Temp);
% hold on
% yline(2.71,'-',{'Chi-Square Threshold'});
% hold off
% xlabel('m_A_i_r  (m/s)');
% ylabel('Negative Log Likelihood Ratio');
% xlim([0 3000]);
% x1 = interp1(LLRatio_Temp, soln_space, 0);
% zero_point = find(soln_space == x1);
% x2 = interp1(LLRatio_Temp(1:zero_point), soln_space(1:zero_point), 2.71);
% x3 = interp1(LLRatio_Temp(zero_point:end), soln_space(zero_point:end), 2.71);
% xline(x1, '--', {'Optimal Parameter Value'});
% xline(x2, '--', {'90% Confidence Interval'});
% xline(x3, '--', {'90% Confidence Interval'});
% ax = gca;
% ax.FontSize = font_size;


%% MAPE
% forecast1 = v.T_PT;
% observed1 = u.T_PT(t);
% ABS1 = abs((observed1 - forecast1)./forecast1)*100;
% MAPE1 = 1/length(t)*sum(ABS1);
% 
% forecast2 = v.L_PT;
% observed2 = u.L_PTA(t);
% ABS2 = abs((forecast2 - observed2)./observed2)*100;
% MAPE2 = 1/length(t)*sum(ABS2);
% 
% MAPEPT = [MAPE1 MAPE2]

%% MPC Initialisation

% function that generates the optimal sequence of
% actions given the currrent starting state and SP
% (initialise and run for first time step)

options = optimoptions('fmincon','Display','off');

sol.y = [p.L_SS; x0.h_PT]; % Ground truth over the first time interval

z.L = []; z.F_in = [];             % Initialise measurement structure
z = Meas(sol, u, 0, p, z); % Load in the desired measurements 
                                   % (choose to use plant data or generated 
                                   % noisy measurements in the Meas function)
x   = z;   % Set the state estimates to initially be equal to the measurements
x.h_PT = x0.h_PT;
x.P = 0.1; % Set the initial covariance

u_opt = fmincon(@(uMV) cost(t(1), uMV, u, p, s, x), p.uvec_init,...
    [], [], [], [],p.MV_min,p.MV_max, [], options); % Everything goes into fmincon,
                                                    % with constraints on
                                                    % the MV movement
                                                    
MV        = u_opt(1); % Set the MV to the first optimal value
output.MV = @(t) MV; 

sol = ode45(@(t,x) PreCoolingTowersODEs(s,p,x,u,t,output,z), [t(1) t(2)], sol.y); % Calculate the ground truth for the first time interval
response(1,:) = deval(sol, t(1));
x.h_PT = sol.y(2,end);

for b = 1:2
    saved.SP_init(b,:) = p.SP; % Save SPs
end

%% MPC Loop

for i = p.loop % Using shorter loop for testing to reduce run-time
        
        % Set the new state to equal the previous prediction
        z = Meas(sol, u, t(i), p, z);
        x = KalmanFilterPT(x, s, output, u, z, [t(i-1) t(i)], p);

        % Change SP and save it
        for j = 1:1:size(p.SP_times,1)
            if p.SP_times(j) == i*p.Ts
                p.SP = p.SP_samples(j)*ones(1,p.N);
            else
                p.SP = p.SP;
            end
        end
        saved.SP_loop(i,:) = p.SP; % Save SPs

        % Perform optimisation
        u_opt = fmincon(@(uMV) cost(t(i), uMV, u, p, s, x), u_opt, [], [], [], [],...
               p.MV_min,p.MV_max, [], options);
        MV(end+1) = u_opt(1); 
        output.MV = griddedInterpolant(t(1:i), MV, 'previous');
        sol = odextend(sol, @(t,x) PreCoolingTowersODEs(s, p, x, u, t, output,z), t(i+1)); % Ground truth (how the system is actually responding)
        response(i,:) = deval(sol, t(i));
      
        fprintf('%d\n',i) 
        x.h_PT(end+1) = sol.y(2,end);
end
v = PTIntermediates(x,u,p,t);
%% MPC Results

saved.SP(1:size(saved.SP_init,1)) = saved.SP_init(:,1);
saved.SP(size(saved.SP_init,1)+1:length(saved.SP_loop)+1) = saved.SP_loop(2:end,1);
saved.SP = saved.SP(1:end-1);

lower_limit = ones(1,length(p.loop)+1).*p.PV_min(1);
upper_limit = ones(1,length(p.loop)+1).*p.PV_max(1);

ax1 = subplot(2,2,1);
hold on
title('PV');
plot(t(1:p.loop(end))/p.Ts, response(:,1),'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, saved.SP', 'r', 'LineWidth',0.5);
hold on
plot(t(1:p.loop(end))/p.Ts, x.L, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.L, '.k', 'MarkerSize', 5);
hold on
plot(t(1:p.loop(end))/p.Ts, lower_limit, '--k', 'MarkerSize', 3);
hold on
plot(t(1:p.loop(end))/p.Ts, upper_limit, '--k', 'MarkerSize', 3);
hold on
ylim([0 110]);
ylabel('Level_P_T (%)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit');
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:i)/p.Ts, output.MV.Values(1:i),'b','LineWidth',1.5);
hold on
ylim([0 900]);
hold off
ylabel('F_o_u_t_P_T (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
title('DV');
plot(t(1:p.loop(end))/p.Ts, u.F_in_generated(t(1:p.loop(end))),'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.F_in, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.F_in, '.k', 'MarkerSize', 5);
hold on
ylim([0 600]);
ylabel('F_i_n_P_T (L/S)'); xlabel('Time (min)');
legend('Ground Truth', 'State Estimate', 'Measurement');
hold off

% ax4 = subplot(2,2,4);
% plot(t(1:p.loop(end))/p.Ts, response(:,2),'color',[0 0.5 0], 'LineWidth',1.5); 
% hold on
% plot(t(1:p.loop(end))/p.Ts, x.h_PT, 'c', 'LineWidth',1);
% ylabel('h_P_T (kJ/kg)'); xlabel('Time (min)');
% legend('Ground Truth', 'State Estimate');
% hold off

ax4 = subplot(2,2,4);
hold on
title('Outlet Temperature');
plot(t(1:p.loop(end))/p.Ts, v.T_PT,'.','Color',[0.8500 0.3250 0.0980]); 
hold on
ylim([12 20]);
ylabel('T_P_T (^oC)'); xlabel('Time (min)');
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');