%% MPC for Level Control: Square Dam
%  Jordan Anastasiou, 10-2023
%  This code is for MPC of level in the
%  square dam, with a KF for state estimation
clc
clear
clf

%% Load Exogenous Inputs
% Loads the measurements as gridded interpolants.
load SavedInterpolantsSD.mat

%% Define process parameters
% Process
p.rho_Water = 1000;    % kg/m3, Density of water
p.m_SDmax   = 6000000; % kg, Maximum mass capacity of SD
p.height_SD = 3;       % m,  Height of the SD
p.area_SD   = 2000;    % m2, Surface area of the SD

p.regressedparameterfields = {'m_evapSD'};      % Field name for the regressed parameter
p.m_evapSD =  1;                                % kg/s, Rate of evaporation from SD (initial guess)
pmEvapVec = S2V(p, p.regressedparameterfields); % Convert the unknown parameter to a vector

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
p.Q = diag([0.09 ; p.w_K(2)]);  % Process noise covariance matrix based on variance from AR model 
% p.A = [1 p.C_L;    
%        0 a.F_in];                 % Transition matrix
p.A = [0 p.C_L;    
       0 0];                        % Transition matrix


%% Define MPC parameters
p.L_SS       = 70;                          % %, Steady state level in the square dam (initial PV value)
p.N          = 5;                           % ~, Number of samples of the control input/prediction nodes
p.Ts         = 60;                          % s, Sampling period, the frequency at which a new control input is determined
p.Stp        = length(t);                   % ~, Number of steps in simulation
p.TL         = (p.Ts*p.Stp) - p.Ts;         % s, Total time or Time Limit (sampling period x total time)
p.loop       = 2:1:200;                     % ~, Simulation points for the loop (mostly for testing)
p.F_outSS    = u.F_outSD(1);                % L/s, Steady state outlet flowrate (initial MV value)
p.uvec_init  = p.F_outSS*ones(1,p.N);       % L/s, Initial points (sequence guess)
p.SP         = 85*ones(1,p.N);              % %, Initial SP for the level in the square dam (initial SP value)
p.SP_changes = 2296;                        % ~, Number of SP changes (2296 - every 22 min. 1335 - every 38 min)
p.SP_min     = 75;                          % %, Lowest SP for the level in the Dam
p.SP_max     = 95;                          % %, Highest SP for the level in the Dam
p.SP_samples = p.SP_min...
               + (p.SP_max - p.SP_min)...
               * rand(p.SP_changes, 1);     % Sample SP changes
p.SP_times   = (0:p.TL/p.SP_changes:p.TL)'; % Times at which the SP should change
p.MV_min     = 0*ones(1,p.N);               % L/s, Minimum MV limit
p.MV_max     = 650*ones(1,p.N);             % L/s, Maximum MV limit
p.PV_min     = 40*ones(1,p.N);              % %, Minimum PV limit (minimum level)
p.PV_max     = 95*ones(1,p.N);              % %, Maximum PV limit (maximum level)
p.Q_Weight   = 0.1;                         % SP weight
p.R_Weight   = 1;                           % MV weight

%% Define state structure and initial conditions
s.statefields    = {'L'};  % Field names for each state 
s.MPCstatefields = {'L','F_in'};
s.KFstatefields  = {'L', 'F_in', 'P'}; % Field names for each state 
x0.L = u.L_SD(0);                      % %, Initial value for level in SD
x0.F_in = u.F_inSD(0);                 % L/s, Initial value for inlet flowrate
x0.P = 0.1; % initial state estimate covariance matrix
x0_vec = S2V(x0, s.statefields);

%% Generate variable inlet flowrate
% In the ARIMA Models file, autoregressive data is produced for use in the
% KF. An additional set is created using the same model to be used in the 
% MPC. This has been named u.F_in_generated.

u.F_in_generated(t);
%plot(t, u.F_in_generated(t));


%% Simulate system of ODEs
% [~, x_vec] = ode45(@(t, x) SquareDamODEs(s, p, x, u, t, output), t, x0_vec);
%x = V2S(x_vec', s.statefields);
% v = SDIntermediates(x, u, p, t);

%% Plot

% Measured vs calculated. Plot the variables for which there is both
% measured and calculated data. These are also the variables used in
% the objective function.
% figure (1)
% plot(t, u.L_SD(t), t, v.L_SD) % %,  Level in SD
% legend('measured', 'predicted');
% xlabel('Time (s)');
% ylabel('L_S_D (%)');

font_size = 17;
% %% Regression
% 
% options = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective',...
%           'CheckGradients', true, 'Display','iter-detailed',...
%           'FiniteDifferenceType','forward', 'StepTolerance', 1e-15,...
%           'FiniteDifferenceStepSize', 0.001, 'MaxFunctionEvaluations', 200,...
%           'MaxIterations',200);
% p_est    = lsqnonlin(@(pmEvapVec) SDCalcError(pmEvapVec, u, p, s, t), pmEvapVec)
% 
% 
% [E, x, v] = SDCalcError(p_est, u, p, s, t);
% 
% % Plot results using parameter estimate/regressed parameter
% figure (2)
% plot(t, u.L_SD(t), t, v.L_SD) % %,  Level in SD
% legend('measured', 'predicted');
% xlabel('Time (s)');
% ylabel('L_S_D (%)');
% 
% %save SquareDam.mat
% 
% %% Likelihood profile
% SSR_Level = []; % Sum of squared residuals for level
% soln_space = 0.001:0.1:4;
% for pEvapvec = soln_space
%     [E, x, v] = SDCalcError(pEvapvec, u, p, s, t);
%     E_2 = E.^2;
%     sum_E_2_Level = sum(E_2(:,1));
%     SSR_Level = [SSR_Level sum_E_2_Level];
% end
% font_size = 17;
% figure(3)
% title('Likelihood Profile for parameter Evap and Level');
% LLRatio_Level = 2*log(SSR_Level/min(SSR_Level));
% plot(soln_space, LLRatio_Level);
% hold on
% yline(2.71,'-',{'Chi-Square Threshold'});
% hold off
% xlabel('m_e_v_a_p (kg/s)');
% ylabel('Negative Log Likelihood Ratio');
% xlim([0 2]);
% x1 = interp1(LLRatio_Level, soln_space, 0);
% zero_point = find(soln_space == x1);
% x2 = interp1(LLRatio_Level(1:zero_point), soln_space(1:zero_point), 2.71);
% x3 = interp1(LLRatio_Level(zero_point:end), soln_space(zero_point:end), 2.71);
% xline(x1, '--', {'Optimal Parameter Value'});
% xline(x2, '--', {'90% Confidence Interval'});
% xline(x3, '--', {'90% Confidence Interval'});
% ax = gca;
% ax.FontSize = font_size;

%% MPC Initialisation

% function that generates the optimal sequence of
% actions given the currrent starting state and SP
% (initialise and run for first time step)

options = optimoptions('fmincon','Display','off');
sol.y = p.L_SS; % Ground truth over the first time interval

z.L = []; z.F_in = [];             % Initialise measurement structure
z = Meas(sol, u, 0, p, z); % Load in the desired measurements 
                                   % (choose to use plant data or generated 
                                   % noisy measurements in the Meas function)
x   = z;   % Set the state estimates to initially be equal to the measurements
x.P = 0.1; % Set the initial covariance

u_opt = fmincon(@(uMV) cost(t(1), uMV, u, p, s, x), p.uvec_init,...
    [], [], [], [],p.MV_min,p.MV_max, [], options); % Everything goes into fmincon,
                                                    % with constraints on
                                                    % the MV movement
                                                    
MV        = u_opt(1); % Set the MV to the first optimal value
output.MV = @(t) MV; 


sol = ode45(@(t,x) SquareDamODEs(s,p,x,u,t,output,z), [t(1) t(2)], sol.y); % Calculate the ground truth for the first time interval
response(1,:) = deval(sol, t(1));

for b = 1:2
    saved.SP_init(b,:) = p.SP; % Save SPs
end

%% MPC Loop

for i = p.loop % Using shorter loop for testing to reduce run-time
        
        % Set the new state to equal the previous prediction
        z = Meas(sol, u, t(i), p, z);
        x = KalmanFilterSD(x, s, output, u, z, [t(i-1) t(i)], p);

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
        sol = odextend(sol, @(t,x) SquareDamODEs(s, p, x, u, t, output,z), t(i+1)); % Ground truth (how the system is actually responding)
        response(i,:) = deval(sol, t(i));
      
        fprintf('%d\n',i)       
end

%% MPC Results

saved.SP(1:size(saved.SP_init,1)) = saved.SP_init(:,1);
saved.SP(size(saved.SP_init,1)+1:length(saved.SP_loop)+1) = saved.SP_loop(2:end,1);
saved.SP = saved.SP(1:end-1);

lower_limit = ones(1,length(p.loop)+1).*p.PV_min(1);
upper_limit = ones(1,length(p.loop)+1).*p.PV_max(1);

ax1 = subplot(3,1,1);
plot(t(1:p.loop(end))/p.Ts, response,'color',[0 0.5 0], 'LineWidth',1.5); 
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
ylim([0 100]);
ylabel('Liquid level (%)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit');
hold off

ax2 = subplot(3,1,2);
plot(t(1:i)/p.Ts, output.MV.Values(1:i),'b','LineWidth',1.5);
hold on
hold off
ylabel('F_o_u_t (L/S)'); xlabel('Time (min)');

ax3 = subplot(3,1,3);
plot(t(1:p.loop(end))/p.Ts, u.F_in_generated(t(1:p.loop(end))),'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.F_in, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.F_in, '.k', 'MarkerSize', 5);
hold on
ylim([0 800]);
ylabel('F_i_n (L/S)'); xlabel('Time (min)');
legend('Ground Truth', 'State Estimate', 'Measurement');
hold off

linkaxes([ax1,ax2,ax3],'x');