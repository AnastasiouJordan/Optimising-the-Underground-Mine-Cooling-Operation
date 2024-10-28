%% Initialise

clc
clear

load SavedInterpolantsCD.mat
load ArimaModelsCD.mat

%% Kalman Filter on Level - KF1
% Taking a balance over the Chill Dam and RPs with the PT outlet, Ice Plant
% flow and UG flow being the incoming and outgoing flows

% States: x = [x1; x2; x3; x4]
% x1 = L_CD
% x2 = F_outPT as the inlet flowrate
% x3 = F_Ice as the first outlet flowrate
% x4 = F_UG as the second outlet flowrate

% Define input data
% This can either be the raw data from the plant, saved and loaded from the
% AR model file as Raw_L', Raw_F_outPT', Raw_F_Ice', Raw_F_UG'
% Or it can be generated data from the AR model file, saved and loaded as
% g_L1, g_F_outPT, g_F_Ice, g_F_UG

% Meas_L       = g_L1;  
% Meas_F_outPT = g_F_outPT;
% Meas_F_Ice   = g_F_Ice;
% Meas_F_UG    = g_F_UG;
Meas_L       = Raw_L';  
Meas_F_outPT = Raw_F_outPT';
Meas_F_Ice   = Raw_F_Ice';
Meas_F_UG    = Raw_F_UG';



xhat_0_1 = [70; 400; 0; 200];       % initial state estimates
P0_1 = eye(length(xhat_0_1))*0.01;  % initial state estimate covariance matrix
H_1 = eye(length(xhat_0_1));        % measurement/observation matrix 


% Define the constant and error matrices
c_K_1 = [0 ; c_F_outPT ; c_F_Ice ; c_F_UG];   % Constants in AR model
w_K_1 = [w_L ; w_F_outPT ; w_F_Ice ; w_F_UG]; % Variance from AR model

% Define the noise
% Measurement noise from instrumentation data sheets
L_sensor_noise_1    = 0.0000015;     % %, Level sensor error
F_outPT_meter_noise = 0.0000001;     % %, Inlet stream flowmeter error
F_Ice_meter_noise   = 0.00001;     % %, Outlet stream flowmeter error
F_UG_meter_noise    = 0.00001;     % %, Outlet stream flowmeter error


L_noise_1     = L_sensor_noise_1   /100*max(Meas_L);        % m,   Std dev of level measurement noise
F_outPT_noise = F_outPT_meter_noise /100*max(Meas_F_outPT); % L/s, Std dev of inlet flowrate measurement noise
F_Ice_noise   = F_Ice_meter_noise/100*max(Meas_F_Ice);      % L/s, Std dev of outlet flowrate measurement noise
F_UG_noise    = F_UG_meter_noise/100*max(Meas_F_UG);        % L/s, Std dev of outlet flowrate measurement noise
R_1           = diag([L_noise_1^2 ;...
                    F_outPT_noise^2 ;...
                    F_Ice_noise^2 ;...
                    F_UG_noise^2]);                 % measurement noise covariance (measurement error squared)
                                                    % Note that the
                                                    % measurement noise has
                                                    % been increased by a
                                                    % factor of 10 in each
                                                    % case to decrease the
                                                    % Kalman gain and to
                                                    % account for older
                                                    % equipment that may
                                                    % not be as accurate.

Var_1 = diag([w_K_1(1) ; w_K_1(2) ; w_K_1(3) ; w_K_1(4)]); % total covariance based on variance from AR model
%Q_1 = (Var_1 - R_1);                                            % process noise covariance 
Q_1 = diag([w_K_1(1) ; w_K_1(2) ; w_K_1(3) ; w_K_1(4)]); % total covariance based on variance from AR model
                                                         % process noise for level = 0 (certain of model)
                          
                                                 
% Define the transition matrix
A_1 = [1 C_L -C_L -C_L;...
       0 a_F_outPT 0 0;...
       0 0 a_F_Ice 0;...
       0 0 0 a_F_UG];


% Define the state estimates vs measurement matrices

x_K_1 = [k_L2 ; k_F_outPT ; k_F_Ice ; k_F_UG]; % State vector estimate/model prediction
z_K_1 = [Meas_L ; Meas_F_outPT ; Meas_F_Ice ; Meas_F_UG]; % Observation vector/measurements

P_1 = P0_1;  % Setting initial state estimate covariance matrix values 


 for n = 2:length(t)-1

   % Priori predictions
   x_K_1(:,n) = A_1*x_K_1(:,n-1) + c_K_1; % priori predicted state
% 
%    if ((x_K_1(3,n) <= 0) || (x_K_1(3,n) >= 800))         % If the predicted Ice Plant flowrate is below zero or above 700 L/s,
%        Q_1 = diag([w_K_1(1) ; w_K_1(2) ; 0 ; w_K_1(4)])      % Set the process noise on the Ice Plant flowrate to zero
%    else
%        Q_1 = diag([w_K_1(1) ; w_K_1(2) ; w_K_1(3) ; w_K_1(4)]) % Otherwise, let it be high so as to track the measurement.
%    end  
   P_1 = A_1.*P_1.*A_1' + Q_1;              % priori predicted covariance

   % Kalman gain
   K_1 = P_1*H_1'*inv(H_1*P_1*H_1' + R_1);
   
   
   % Correction
   e_1(:,n) = z_K_1(:,n) - (H_1*x_K_1(:,n)); % Measurement residual where z is the observation
                                     % vector or measurements
                               
   x_K_1(:,n) = x_K_1(:,n) + (K_1*e_1(:,n)); % posteriori state estimate
   P_1 = P_1 - (K_1*H_1*P_1);                  % posteriori covariance

end

x_K_L1      = x_K_1(1,:); % Kalman estimate for Chill Dam Level
x_K_F_outPT = x_K_1(2,:); % Kalman estimate for Chill Dam Inlet Flowrate
x_K_F_Ice   = x_K_1(3,:); % Kalman estimate for Chill Dam Outlet Flowrate
x_K_F_UG    = x_K_1(4,:); % Kalman estimate for Chill Dam Outlet Flowrate


% Plot the results of the Kalman Estimate, AR Model Estimate and the Raw
% Data
figure(1)
title('Kalman Filter');
subplot(4,1,1)
plot(t, x_K_L1', t, Meas_L, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('L_C_D (m)')
subplot(4,1,2)
plot(t, x_K_F_outPT', t, Meas_F_outPT, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('F_o_u_t_P_T (L/s)')
subplot(4,1,3)
plot(t, x_K_F_Ice', t, Meas_F_Ice, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('F_I_c_e (L/s)')
subplot(4,1,4)
plot(t, x_K_F_UG', t, Meas_F_UG, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('F_U_G (L/s)')

u.L_CD1_filtered    = griddedInterpolant(t, x_K_L1');
u.F_out_PT_filtered = griddedInterpolant(t, x_K_F_outPT');
u.F_Ice_filtered    = griddedInterpolant(t, x_K_F_Ice');
u.F_UG_filtered     = griddedInterpolant(t, x_K_F_UG');

%% Kalman Filter on Level - KF2
% Taking a balance over the Chill Dam alone
% States: x = [x1; x2; x3]
% x1 = L
% x2 = Fin
% x3 = Fout

% Define input data
% This can either be the raw data from the plant, saved and loaded from the
% AR model file as Raw_L', Raw_F_in', Raw_F_out'
% Or it can be generated data from the AR model file, saved and loaded as
% g_L2, g_F_in, g_F_out

% Meas_L       = g_L2;  
% Meas_F_in    = g_F_in;
% Meas_F_out   = g_F_out;
Meas_L       = Raw_L';  
Meas_F_in    = Raw_F_in';
Meas_F_out   = Raw_F_out';


xhat_0_2 = [70; 400; 350];        % initial state estimates
P0_2 = eye(length(xhat_0_2))*0.01;  % initial state estimate covariance matrix
H_2 = eye(length(xhat_0_2));        % measurement/observation matrix 


% Define the constant and error matrices
c_K_2 = [0 ; c_F_in ; c_F_out];   % Constants in AR model
w_K_2 = [w_L ; w_F_in ; w_F_out]; % Variance from AR model

% Define the noise
% Measurement noise from instrumentation data sheets
L_sensor_noise_2  = 0.05;       % %, Level sensor error
F_in_meter_noise  = 0.009;        % %, Inlet stream flowmeter error
F_out_meter_noise = 0.031;        % %, Outlet stream flowmeter error


L_noise_2   = L_sensor_noise_2   /100*max(Meas_L);   % m,   Std dev of level measurement noise
F_in_noise  = F_in_meter_noise /100*max(Meas_F_in);  % L/s, Std dev of inlet flowrate measurement noise
F_out_noise = F_out_meter_noise/100*max(Meas_F_out); % L/s, Std dev of outlet flowrate measurement noise
R_2         = diag([L_noise_2^2 ;...
                    F_in_noise^2 ;...
                    F_out_noise^2]);                % measurement noise covariance (measurement error squared)
                                                    % Note that the
                                                    % measurement noise has
                                                    % been increased by a
                                                    % factor of 10 in each
                                                    % case to decrease the
                                                    % Kalman gain and to
                                                    % account for older
                                                    % equipment that may
                                                    % not be as accurate.

Var_2 = diag([w_K_2(1) ; w_K_2(2) ; w_K_2(3)]); % total covariance based on variance from AR model
%Q_2 = (Var_2 - R_2);                            % process noise covariance 
Q_2 = diag([w_K_2(1) ; w_K_2(2) ; w_K_2(3)]); % total covariance based on variance from AR model
                                              % process noise for level = 0 (certain of model)
                          
                                                 
% Define the transition matrix
A_2 = [1 C_L -C_L;...
       0 a_F_in 0;...
       0 0 a_F_out];


% Define the state estimates vs measurement matrices
x_K_2 = [k_L1 ; k_F_in ; k_F_out]; % State vector estimate/model prediction
z_K_2 = [Meas_L ; Meas_F_in ; Meas_F_out]; % Observation vector/measurements

P_2 = P0_2;  % Setting initial state estimate covariance matrix values 

e_2 = [];
for n = 2:length(t)-1
   % Priori predictions
   x_K_2(:,n) = A_2*x_K_2(:,n-1) + c_K_2; % priori predicted state
   P_2 = A_2.*P_2.*A_2' + Q_2;              % priori predicted covariance

   % Kalman gain
   K_2 = P_2*H_2'*inv(H_2*P_2*H_2' + R_2);
   

   % Correction
   e_2(:,n) = z_K_2(:,n) - (H_2*x_K_2(:,n)); % Measurement residual where z is the observation
                                             % vector or measurements
                               
   x_K_2(:,n) = x_K_2(:,n) + (K_2*e_2(:,n)); % posteriori state estimate
   P_2 = P_2 - (K_2*H_2*P_2);                  % posteriori covariance
end

x_K_L2    = x_K_2(1,:); % Kalman estimate for Chill Dam Level
x_K_F_in  = x_K_2(2,:); % Kalman estimate for Chill Dam Inlet Flowrate
x_K_F_out = x_K_2(3,:); % Kalman estimate for Chill Dam Outlet Flowrate

% Plot the results of the Kalman Estimate, AR Model Estimate and the Raw
% Data
figure(2)
title('Kalman Filter');
subplot(3,1,1)
plot(t, x_K_L2', t, Meas_L, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('L_C_D (m)')
subplot(3,1,2)
plot(t, x_K_F_in', t, Meas_F_in, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('F_i_n_C_D (L/s)')
subplot(3,1,3)
plot(t, x_K_F_out', t, Meas_F_out, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('F_o_u_t_C_D (L/s)')

u.L_CD2_filtered = griddedInterpolant(t, x_K_L2');
u.F_in_filtered  = griddedInterpolant(t, x_K_F_in');
u.F_out_filtered = griddedInterpolant(t, x_K_F_out');

save KalmanFilterCD.mat u
