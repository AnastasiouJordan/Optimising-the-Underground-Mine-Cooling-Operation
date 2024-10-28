%% System of ODEs: Square Dam
%  Jordan Anastasiou, 07-2022
%  This code is for the square dam
clc
clear
clf

%% Define parameters
p.rho_Water = 1000;       % kg/m3, Density of water
p.m_SDmax   = 6000000;    % kg,    Maximum mass capacity of SD

p.regressedparameterfields = {'m_evapSD'};

% Initial guess for evporation rate
p.m_evapSD =  1; % kg/s, Rate of evaporation from SD

pmEvapVec = S2V(p, p.regressedparameterfields); % Convert the unknown parameter to a vector

%% Define exogenous inputs

load SavedInterpolantsSD.mat

%% Define state structure and initial conditions
s.statefields = {'m_SD'};          % Field names for each state 
x0.m_SD = u.L_SD(0)*p.m_SDmax/100; % kg, Initial value for mass in SD
x0_vec = S2V(x0, s.statefields);

%% Load Kalman Filtered Values
load KalmanFilterSD.mat

%% Define MPC parameters
control.L_SS = x0_vec; % kg, Steady state mass in the square dam


% Define the prediction horizon
control.T = 10; % Final time of prediction horizon
control.N = 10; % Number of samples of the control input

% Define the sampling time
control.Ts  = 60;  % s, Sampling period: the frequency at which a new control input is determined
control.Stp = length(t); % Number of time steps to be simulated (length of simulation)
control.TL  = (control.Ts*control.Stp) - control.Ts; % Total time or Time Limit (sampling period x total time)
control.HS  = 1; % Horizon start (first time step in prediction horizon)

% Define the steady state values
control.F_inSS = u.F_in_filtered.Values(1); % Define the DV from the model
control.uvec_init = control.F_inSS*ones(1,control.N); % initial points (sequence guess)

% Define the SP and DVs
control.SP         = 80*p.m_SDmax/100; % kg, SP for the mass in the square dam
control.SP_changes = 25543; % number of SP changes
control.SP_min     = 40*p.m_SDmax/100; % kg, Lowest SP for the mass in the Dam
control.SP_max     = 80*p.m_SDmax/100; % kg, Highest SP for the mass in the Dam
control.SP_samples = control.SP_min...
                     + (control.SP_max - control.SP_min)...
                     * rand(control.SP_changes, 1); % sample SP changes
control.SP_times   = (0:control.TL/control.SP_changes:control.TL)'; % Times at which the SP should change


control.F_in       = u.F_in_filtered(t); % Steady state inlet flow rate (L/s)
control.DV_changes = 0; % Number of DV changes (changes in inlet flowrate)
control.DV_min     = 200; % L/s, Lowest inlet flowrate
control.DV_max     = 1000; % L/s, Highest inlet flowrate
control.DV_samples = control.DV_min...
                     + (control.DV_max - control.DV_min)...
                     * rand(control.DV_changes,1); % sample DV changes
control.DV_times   = randi(control.Stp, control.DV_changes, 1);

control.gamma = 0.1; % discount factor used in the discounted MPC cost function, 2023-03-30
control.R_bound_high = -1*-0.25; % lower bound for squared error scaling (2023-04-01)
control.R_bound_low = -1*0; % upper bound for squared error scaling (2023-04-01)

control.MPC_ref_undisc = 0; % initialize undiscounted cost achieved on the "real process" by the MPC (2023-04-01)
control.MPC_ref_disc = 0; % initialize discounted cost achieved on the "real process" by the MPC (2023-04-01)

%% Simulate system of ODEs
[~, x_vec] = ode45(@(t, x) SquareDamODEs(s, p, x, u, t), t, x0_vec);
x = V2S(x_vec', s.statefields);
v = SDIntermediates(x, u, p, t);

%% Plot
font_size = 18;
%Measured vs calculated. Plot the variables for which there is both
%measured and calculated data. These are also the variables used in
%the objective function.
figure (1)
plot(t/86400, u.L_SD(t), t/86400, v.L_SD) % %,  Level in SD
legend('measured', 'predicted');
ylabel('L_S_D (%)');
xlabel('Time (days)');
ax = gca;
ax.FontSize = font_size;

%% Regression

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

%save SquareDam.mat


%% Controller

% 
% %for i = 1:1:control.TL/control.Ts
% for i = 1:1:200
%     % change SP if necessary
%     for j = 1:1:size(control.SP_times,1)
%         if control.SP_times(j) == i
%             control.SP = control.SP_samples(j);
%         end
%     end
% %     % change DV if necessary
% %     for j = 1:1:size(control.DV_times,1)
% %         if control.DV_times(j) == i
% %             control.F_in = control.DV_samples(j);
% %         end
% %     end
% 
%     if i == 1
%         control.L_init = control.L_SS; % Set the initial level to the steady-state level specified
%         [u_opt,model,t] = optimise_traj(t, u, p, control); % obtain optimal sequence of control inputs across the control horison
%         control.model = u_opt(1); % select initial control input from optimized trajectory
%         %p.x0 = p.x; % set starting valve position for next optimization step
%         TSPAN = i:60:i+control.Ts; %i:1:i+p.sampling_period - 1; (2023-03-31)
%         [time, x_output] = ode45(@(time, x) SquareDamODEs(s, p, x, u, time), TSPAN, control.L_SS);
%         prevTime = time(end); % final time (min)
%         prevH = x_output(end); % final liquid height (m)
%         saved_trajectories.T(:,i+1) = time(end); % save time trajectory (2023-03-31)
%         saved_trajectories.L(:,i+1) = x_output(end); % save liquid height trajectory (2023-03-31)
%         saved_trajectories.SP(:,i) = control.SP; % save SPs
%         saved_trajectories.F_in(:,i) = control.F_in'; % save DVs
%         saved_trajectories.MV(:,i) = control.model; % control input
%         control.MPC_ref_undisc = control.MPC_ref_undisc + ( ( (control.SP - x_output(end))^2 - control.R_bound_low )/( control.R_bound_high - control.R_bound_low ) ); % update MPC cost on the "real process" (2023-04-01)
%         control.MPC_ref_disc = control.MPC_ref_disc + ( ( (control.SP - x_output(end))^2 - control.R_bound_low )/( control.R_bound_high - control.R_bound_low ) ); % update MPC cost on the "real process" (2023-04-01)
% 
%     elseif i > 1
%         control.L_init = prevH;%crntstartH; % set starting valve position for next optimization step (2023-03-31)
%         [u_opt, model] = optimise_traj(t, u, p, control); % obtain optimal sequence of control inputs across the control horison
%         control.model = u_opt(1); % select initial control input from optimized trajectory
%         TSPAN = i:60:i+control.Ts;%prevTime+1:1:prevTime + p.sampling_period;%prevTime:1:prevTime+p.sampling_period - 1; (2023-03-31)
%         %prevH = HOutput(end); % final liquid height (m)
%         [time, x_output] = ode45(@(time, x) SquareDamODEs(s, p, x, u, time), TSPAN, prevH);
%         prevTime = time(end); % final time (min)
%         prevH = x_output(end); % final liquid height (m)
%         saved_trajectories.T(:,i+1) = time(end); % save time trajectory (2023-03-31)
%         saved_trajectories.H(:,i+1) = x_output(end); % save liquid height trajectory (2023-03-31)
%         saved_trajectories.SP(:,i) = control.SP; % save SPs
%         saved_trajectories.F_in(:,i) = control.F_in'; % save DVs
%         saved_trajectories.MV(:,i) = control.model; % control input
%         control.MPC_ref_undisc = control.MPC_ref_undisc + ( 1^(i-1) )*( ( (control.SP - x_output(end))^2 - control.R_bound_low )/( control.R_bound_high - control.R_bound_low ) ); % update MPC cost on the "real process" (2023-04-01)
%         control.MPC_ref_disc = control.MPC_ref_disc + ( control.gamma^(i-1) )*( ( (control.SP - x_output(end))^2 - control.R_bound_low )/( control.R_bound_high - control.R_bound_low ) ); % update MPC cost on the "real process" (2023-04-01)
% 
%     end
% 
%     fprintf('%d\n',i)
% end
% 
% figure(2)
% saved_trajectories.SP = saved_trajectories.SP/p.m_SDmax*100; % Converting back to level as a percentage
% saved_trajectories.H  = saved_trajectories.H/p.m_SDmax*100; % Converting back to level as a percentage
% 
% plot(saved_trajectories.SP,'k:','LineWidth',2); hold on; 
% plot(saved_trajectories.H,'k--','LineWidth',2); set(gcf,'color','white'); hold on 
% xlabel('Time (s)'); ylabel('Level (%)'); hold on
% set(gca,'fontsize',20)
% legend('SP','NMPC'); axis tight;
% 
% 
% figure(3)
% subplot(2,1,1)
% plot(linspace(0,control.T,size(model,1)),model(:,1),'r--','LineWidth',3); hold on
% %plot(linspace(0,T,size(x,1)),x(:,2),'b--','LineWidth',3); hold on
% xlabel('Time step in trajectory'); ylabel('model');
% set(gca,'fontsize',20)
% yline(saved_trajectories.SP,'k:','LineWidth',1);
% legend('L (%)','con 1')
% subplot(2,1,2)
% stairs(0:control.T/control.N:control.T,[u_opt,u_opt(end)],'b-','LineWidth',3); hold on;
% xlabel('Time step in trajectory'); ylabel('u'); set(gca,'fontsize',20)
% set(gcf,'color','white');
