%% System of ODEs: Chill Dam
%  Jordan Anastasiou, 2023-11
%  This code is for the chill dam over which the
%  mass and energy balances were performed. The 
%  MPC is built-in for level control, as well as 
%  RP inlet temperature control with the CD outlet 
%  valve as the MV.

clc
clear
clf

%% Define exogeneous inputs
load SavedInterpolantsCD.mat

%% Define Process Parameters
p.regressedparameterfields = {'UA_CD'};

p.rho_Water = 1000;    % kg/m3,  Density of water
p.h_0       = 0.10186; % kJ/kg,  Reference specific enthalpy
p.T_0       = 0.01;    % oC,     Reference temperature
p.C_p       = 4.1831;  % kJ/kgK, Heat capacity of water
p.m_CDmax   = 8000000; % kg,     Maximum mass capacity of CD


p.height_CD = 6.7; % m, Height of the Chill Dam
p.length_CD = 48;  % m, Length of the Chill Dam
p.width_CD  = 25;  % m, Width of the Chill Dam

p.area_CD   = p.length_CD*p.width_CD; % m2, Area of the Chill Dam

% Initial guess for UA_CD

p.UA_CD     = 190; % kJ/Ks,  Heat transfer coefficient for the rate
                   % of heat transfer between the water in the CD and
                  % the environment through the dam walls

pUAvec = S2V(p, p.regressedparameterfields); % Convert the unknown parameter to a vector

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
p.Q = diag([0.009; p.w_K(2)]);  % Process noise covariance matrix based on variance from AR model 
% p.A = [1 p.C_L;    
%        0 a.F_in];                 % Transition matrix
p.A = [0 p.C_L;    
       0 0];                        % Transition matrix

%% Define MPC parameters
p.L_SS       = 80;                          % %, Steady state level in the square dam (initial PV value)
p.N          = 3;                           % ~, Number of samples of the control input/prediction nodes
p.Ts         = 60;                          % s, Sampling period, the frequency at which a new control input is determined
p.Stp        = length(t);                   % ~, Number of steps in simulation
p.TL         = (p.Ts*p.Stp) - p.Ts;         % s, Total time or Time Limit (sampling period x total time)
p.loop       = 2:1:300;                      % ~, Simulation points for the loop (mostly for testing)
p.F_outSS    = u.F_outCD(1);                % L/s, Steady state outlet flowrate (initial MV value)
p.uvec_init  = p.F_outSS*ones(1,p.N);       % L/s, Initial points (sequence guess)
p.SP         = 85*ones(1,p.N);              % %, Initial SP for the level in the square dam (initial SP value)
p.SP_changes = 2296;                        % ~, Number of SP changes (2296 - every 22 min. 1335 - every 38 min)
p.SP_min     = 75;                          % %, Lowest SP for the level in the Dam
p.SP_max     = 90;                          % %, Highest SP for the level in the Dam
p.SP_samples = p.SP_min...
               + (p.SP_max - p.SP_min)...
               * rand(p.SP_changes, 1);     % Sample SP changes
p.SP_times   = (0:p.TL/p.SP_changes:p.TL)'; % Times at which the SP should change
p.MV_min     = 0*ones(1,p.N);               % L/s, Minimum MV limit (ensure that there is always enough flow for RPs)
p.MV_max     = 600*ones(1,p.N);             % L/s, Maximum MV limit
p.PV_min     = 40*ones(1,p.N);              % %, Minimum PV limit (minimum level)
p.PV_max     = 95*ones(1,p.N);              % %, Maximum PV limit (maximum level)
p.Q_Weight   = 1;                           % SP weight
p.R_Weight   = 0.1;                         % MV weight
%% Define state structure and initial conditions

s.statefields    = {'L','h_CD'};              % Field names for each state  
s.MPCstatefields = {'L','F_in', 'h_CD'};      % Field names for each state in MPC
s.KFstatefields  = {'L','F_in', 'h_CD', 'P'}; % Field names for each state in KF

x0.L    = u.L_CD(0);   % %, Initial value level
x0.h_CD = 18;          % kJ/kgK, Initial value for enthalpy of water leaving CD
x0.F_in = u.F_inRP(0)...
        - u.F_UG(0);   % L/s, Initial value for inlet flowrate
x0.P    = 0.1;         % Initial state estimate covariance matrix


x0_vec = S2V(x0, s.statefields);

%% Generate variable inlet flowrate
% In the ARIMA Models file, autoregressive data is produced for use in the
% KF. An additional set is created using the same model to be used in the 
% MPC. This has been named u.F_in_generated.

u.F_in_generated(t);


%% Simulate system of ODEs

% [~, x_vec] = ode45(@(t, x) ChillDamODEs(s, p, x, u, t), t, x0_vec);
% x = V2S(x_vec', s.statefields);
% v = CDIntermediates(x, u, p, t);

%Simulate using ODExtend
% tic
% x_vec = S2V(x0, s.statefields);
% for i = 2:length(t)
%     xk = x_vec(:, i-1);          % Initial values as a vector
%    
%     % Using odextend
%     if i == 2
%         sol = ode45(@(t, x) ChillDamODEs(s, p, x, u, t), [t(i-1) t(i)], xk);
%     else
%         sol = odextend(sol, @(t,x) ChillDamODEs(s, p, x, u, t), t(i), xk);
%     end
%     xnew = sol.y(:, end);
%  
%     x_vec(:, i) = xnew;
% end
% xODEexd = V2S(x_vec, s.statefields);
% toc
% disp('Simulation using odextend complete')
% 
% % Simulate using RK4-steps
% tic
% x_vec = S2V(x0, s.statefields);
% for i = 2:length(t)
%     xk = x_vec(:, i-1);          % Initial values as a vector
%    
%    
%     % Using an RK4 step
%     k1 = ChillDamODEs(s, p, xk,               u, t(i-1));
%     k2 = ChillDamODEs(s, p, xk + deltaT/2*k1, u, t(i-1) + deltaT/2);
%     k3 = ChillDamODEs(s, p, xk + deltaT/2*k2, u, t(i-1) + deltaT/2);
%     k4 = ChillDamODEs(s, p, xk + deltaT*k3,   u, t(i));
%     xnew = (xk + deltaT/6*(k1 + 2*k2 + 2*k3 + k4))';    % Transpose to generate row vector, as would be done using ode45
%  
%     x_vec(:, i) = xnew;
% end
% xRK4 = V2S(x_vec, s.statefields);
% toc
% disp('Simulation using RK4 complete')
% 
% figure(6)
% plot( t, x.h_CD, t, xODEexd.h_CD, t, xRK4.h_CD)
% legend('ODE45','ODExtend', 'RK4')

%% Plot
% Plot results using initial assumption for the parameter value
font_size = 25;
% figure (1)
% title('Prediction Results with Assumed Parameter Values');
% subplot(2,1,1)
% plot(t/86400, u.T_CD(t), t/86400, v.T_CD);  
% legend('measured', 'predicted','FontSize',font_size);
% %text(3100000, 2, sprintf('UA = %.2f kJ/Ks', p.UA_CD));
% xlabel('Time (days)');
% ylabel('T_C_D (^oC)');
% ax = gca;
% ax.FontSize = font_size;
% subplot(2,1,2)
% plot(t/86400, u.L_CD(t), t/86400, v.L_CD);
% legend('measured', 'predicted','FontSize',font_size);
% xlabel('Time (days)');
% ylabel('L_C_D (%)');
% ax = gca;
% ax.FontSize = font_size;

%% Regression

% options = optimoptions('lsqnonlin', 'StepTolerance', 1e-9,...
%                        'Algorithm','trust-region-reflective',...
%                        'FiniteDifferenceType','central',...
%                        'TypicalX', 200,...
%                        'Display', 'iter-detailed');
% 
% p_est    = lsqnonlin(@(pUAvec) CDCalcError(pUAvec, u, p, s, t), pUAvec, 50, 250, options)
% 
% 
% [E, x, v] = CDCalcError(p_est, u, p, s, t);
% 
% 
% % Plot results using parameter estimate/regressed parameter
% 
% figure (2)
% title('Prediction Results with Regressed Parameter Values');
% subplot(2,1,1)
% plot(t, u.T_CD(t), t, v.T_CD);  
% legend('measured', 'predicted');
% text(3100000, 2, sprintf('UA = %.2f kJ/Ks', p.UA_CD));
% xlabel('Time (s)');
% ylabel('T_C_D (^oC)');
% subplot(2,1,2)
% plot(t, u.L_CD(t), t, v.L_CD);
% legend('measured', 'predicted');
% text(3100000, 2, sprintf('UA = %.2f kJ/Ks', p.UA_CD));
% xlabel('Time (s)');
% ylabel('L_C_D (%)');
% ax = gca;
% ax.FontSize = font_size;

% %% Likelihood profile
% SSR_Level = []; % Sum of squared residuals for level
% SSR_Temp  = []; % Sum of squared residuals for temperature
% soln_space = 0.01:10:1000;
% for pUAvec = soln_space
%     [E, x, v] = CDCalcError(pUAvec, u, p, s, t);
%     E_2 = E.^2;
%     sum_E_2_Level = sum(E_2(:,1));
%     sum_E_2_Temp  = sum(E_2(:,2));
%     SSR_Level = [SSR_Level sum_E_2_Level];
%     SSR_Temp  = [SSR_Temp sum_E_2_Temp];
% end
% 
% figure(3)
% title('Likelihood Profile for parameter UA and Level');
% LLRatio_Level = 2*log(SSR_Level/min(SSR_Level));
% plot(soln_space, LLRatio_Level);
% hold on
% yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
% hold off
% xlabel('UA_C_D (kJ/Ks)');
% ylabel('Negative Log Likelihood Ratio');
% xlim([0 1000]);
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
% xlabel('UA_C_D (kJ/Ks)');
% ylabel('Negative Log Likelihood Ratio');
% xlim([0 1000]);
% x1 = interp1(LLRatio_Temp, soln_space, 0);
% zero_point = find(soln_space == x1);
% x2 = interp1(LLRatio_Temp(1:zero_point), soln_space(1:zero_point), 2.71);
% x3 = interp1(LLRatio_Temp(zero_point:end), soln_space(zero_point:end), 2.71);
% xline(x1, '--', {'Optimal Parameter Value'});
% xline(x2, '--', {'90% Confidence Interval'});
% xline(x3, '--', {'90% Confidence Interval'});
% ax = gca;
% ax.FontSize = font_size;
% 
% %% MAPE for all units
% % CD
% forecastTCD = v.T_CD;
% observedTCD = u.T_CD(t);
% ABSTCD = abs((forecastTCD - observedTCD)./observedTCD)*100;
% MAPETCD = 1/length(t)*sum(ABSTCD);
% 
% forecastLCD = v.L_CD;
% observedLCD = u.L_CD(t);
% ABSLCD = abs((forecastLCD - observedLCD)./observedLCD)*100;
% MAPELCD = 1/length(t)*sum(ABSLCD);
% 
% MAPECD = [MAPETCD MAPELCD];
% 
% %% PT
% load PreCoolingTowers.mat
% forecastTPT = v.T_PT;
% observedTPT = u.T_PT(t);
% ABSTPT = abs((observedTPT - forecastTPT)./forecastTPT)*100;
% MAPETPT = 1/length(t)*sum(ABSTPT);
% 
% forecastLPT = v.L_PT;
% observedLPT = u.L_PTA(t);
% ABSLPT = abs((forecastLPT - observedLPT)./observedLPT)*100;
% MAPELPT = 1/length(t)*sum(ABSLPT);
% 
% MAPEPT = [MAPETPT MAPELPT];
% 
% %% RPs
% load FridgePlants.mat
% N = 400:650;
% % RP1
% forecast1 = v.T_RP(:,1);
% observed1 = u.T_RP.Values(N,1);
% ABS1 = abs((forecast1 - observed1)./observed1)*100;
% MAPE1 = 1/length(N)*sum(ABS1);
% 
% % RP2
% forecast2 = v.T_RP(:,2);
% observed2 = u.T_RP.Values(N,2);
% ABS2 = abs((forecast2 - observed2)./observed2)*100;
% MAPE2 = 1/length(N)*sum(ABS2);
% 
% % RP3
% forecast3 = v.T_RP(:,3);
% observed3 = u.T_RP.Values(N,3);
% ABS3 = abs((forecast3 - observed3)./observed3)*100;
% MAPE3 = 1/length(N)*sum(ABS3);
% 
% % RP4
% forecast4 = v.T_RP(:,4);
% observed4 = u.T_RP.Values(N,4);
% ABS4 = abs((forecast4 - observed4)./observed4)*100;
% MAPE4 = 1/length(N)*sum(ABS4);
% 
% % RP5
% forecast5 = v.T_RP(:,5);
% observed5 = u.T_RP.Values(N,5);
% ABS5 = abs((forecast5 - observed5)./observed5)*100;
% MAPE5 = 1/length(N)*sum(ABS5);
% 
% MAPERP = [MAPE1 MAPE2 MAPE3 MAPE4 MAPE5];
% 
% %% SD
% load SquareDam.mat
% forecastLSD = v.L_SD;
% observedLSD = u.L_SD(t);
% ABSLSD = abs((observedLSD - forecastLSD)./forecastLSD)*100;
% MAPELSD = 1/length(t)*sum(ABSLSD);
% 
% 
% %% Plot
% ALL_ABS = [ABSLSD; ABSLPT; ABSTPT; ABS1; ABS2; ABS3; ABS4; ABS5; ABSLCD; ABSTCD];
% group = [repmat({'SD Level'},51087,1); repmat({'PT Level'},51087,1); repmat({'PT Temp'},51087,1); repmat({'RP1 Temp'},251,1);...
%     repmat({'RP2 Temp'},251,1);repmat({'RP3 Temp'},251,1);repmat({'RP4 Temp'},251,1);repmat({'RP5 Temp'},251,1);repmat({'CD Level'},51087,1);repmat({'CD Temp'},51087,1)];
% bh=boxplot(ALL_ABS,group,'PlotStyle','compact','MedianStyle','line','BoxStyle','outline','Symbol','.r','LabelOrientation','horizontal');
% % set(gca,'TickLabelInterpreter','tex');
% % set(gca, 'xticklabel', {'L_S_D', 'L_P_T', 'T_P_T', 'T_R_P_1', 'T_R_P_2','T_R_P_3','T_R_P_4', 'T_R_P_5', 'L_C_D', 'T_C_D'})
% ylim([-2 85]);
% ylabel('Absolute Percentage Error');
% ax = gca;
% ax.FontSize = 30;
% set(bh,'LineWidth', 2);

%% MPC Initialisation

% function that generates the optimal sequence of
% actions given the currrent starting state and SP
% (initialise and run for first time step)

options = optimoptions('fmincon','Display','off');

sol.y = [p.L_SS; x0.h_CD]; % Ground truth over the first time interval

z.L = []; z.F_in = [];             % Initialise measurement structure
z = Meas(sol, u, 0, p, z); % Load in the desired measurements 
                                   % (choose to use plant data or generated 
                                   % noisy measurements in the Meas function)
x   = z;   % Set the state estimates to initially be equal to the measurements
x.h_CD = x0.h_CD;
x.P = 0.1; % Set the initial covariance

u_opt = fmincon(@(uMV) cost(t(1), uMV, u, p, s, x), p.uvec_init,...
    [], [], [], [],p.MV_min,p.MV_max, [], options); % Everything goes into fmincon,
                                                    % with constraints on
                                                    % the MV movement
                                                    
MV        = u_opt(1); % Set the MV to the first optimal value
output.MV = @(t) MV; 

sol = ode45(@(t,x) ChillDamODEs(s,p,x,u,t,output), [t(1) t(2)], sol.y); % Calculate the ground truth for the first time interval
response(1,:) = deval(sol, t(1));
x.h_CD = sol.y(2,end);

for b = 1:2
    saved.SP_init(b,:) = p.SP; % Save SPs
end

%% MPC Loop

for i = p.loop % Using shorter loop for testing to reduce run-time
        
        % Set the new state to equal the previous prediction
        z = Meas(sol, u, t(i), p, z);
        x = KalmanFilterCD(x, s, output, u, z, [t(i-1) t(i)], p);

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
        sol = odextend(sol, @(t,x) ChillDamODEs(s, p, x, u, t, output), t(i+1)); % Ground truth (how the system is actually responding)
        response(i,:) = deval(sol, t(i));
      
        fprintf('%d\n',i) 
        x.h_CD(end+1) = sol.y(2,end);
end
v = CDIntermediates(x,u,p,t(1:p.loop(end)),output);
%% MPC Results

saved.SP(1:size(saved.SP_init,1)) = saved.SP_init(:,1);
saved.SP(size(saved.SP_init,1)+1:length(saved.SP_loop)+1) = saved.SP_loop(2:end,1);
saved.SP = saved.SP(1:end-1);

% Level Limits
lower_limit = ones(1,length(p.loop)+1).*p.PV_min(1);
upper_limit = ones(1,length(p.loop)+1).*p.PV_max(1);

% Fridge plant inlet temperature SP
T_inRP_SP = ones(1,length(p.loop)+1).*12;

font_size = 35;
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
ylabel('Level_C_D (%)'); xlabel('Time (min)');
legend('Ground Truth','SP', 'State Estimate', 'Measurement', 'Lower Limit', 'Upper Limit','Location', 'Best');
ax = gca;
ax.FontSize = font_size;
hold off

ax2 = subplot(2,2,2);
hold on
title('MV');
plot(t(1:i)/p.Ts, output.MV.Values(1:i),'b','LineWidth',1.5);
hold on
ylim([0 700]);
hold off
ax = gca;
ax.FontSize = font_size;
ylabel('F_o_u_t_C_D (L/S)'); xlabel('Time (min)');

ax3 = subplot(2,2,3);
hold on
title('DVs');
colororder([0 0.5 0; 0 0.4470 0.7410])
yyaxis left
plot(t(1:p.loop(end))/p.Ts, u.F_in_generated(t(1:p.loop(end))),'color',[0 0.5 0], 'LineWidth',1.5); 
hold on
plot(t(1:p.loop(end))/p.Ts, x.F_in, 'c', 'LineWidth',1);
hold on
plot(t(1:p.loop(end))/p.Ts, z.F_in, '.k', 'MarkerSize', 5);
hold on
ylim([0 1000]);
ylabel('F_i_n_C_D (L/S)'); xlabel('Time (min)');
hold on
yyaxis right
plot(t(1:p.loop(end))/p.Ts, u.F_Ice_generated((t(1:p.loop(end)))), 'color', [0 0.4470 0.7410], 'LineWidth',1.5)
ylabel('F_I_c_e (L/S)'); xlabel('Time (min)');
legend('Ground Truth', 'State Estimate', 'Measurement','Location', 'Best');
ax = gca;
ax.FontSize = font_size;
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
title('Temperatures');
plot(t(1:p.loop(end))/p.Ts, v.T_CD,'.','Color',[0.8500 0.3250 0.0980]); 
hold on
ylim([0 15]);
plot(t(1:p.loop(end))/p.Ts, v.T_inRP,'.','Color',[0.6350 0.0780 0.1840]);
hold on
plot(t(1:p.loop(end))/p.Ts, T_inRP_SP,'r', 'LineWidth',0.5); 
ylabel('T (^oC)'); xlabel('Time (min)');
legend('T_o_u_t_C_D', 'T_i_n_R_P','T_i_n_R_P SP', 'Location', 'Best');
ax = gca;
ax.FontSize = font_size;
hold off

linkaxes([ax1,ax2,ax3,ax4],'x');
