clc
clear

load ArimaModelsSD.mat
load SavedInterpolantsSD.mat

% CONSTRAINED KALMAN FILTER

% Define the inital conditions
% States: x = [x1; x2; x3]
% x1 = L
% x2 = Fin
% x3 = Fout

% Define input data
% This can either be the raw data from the plant, saved and loaded from the
% AR model file as Raw_L', Raw_F_in', Raw_F_out'.
% Or it can be generated data from the AR model file, saved and loaded as
% g_L, g_F_in, g_F_out

% Comment out the first 3 lines to use raw plant data
% Comment out the last 3 lines to use generated/simulated data

% Meas_L     = g_L;  
% Meas_F_in  = g_F_in;
% Meas_F_out = g_F_out;
Meas_L     = Raw_L';  
Meas_F_in  = Raw_F_in';
Meas_F_out = Raw_F_out';

xhat_0 = [2.9; 400; 400];      % intial state estimates
P0 = eye(length(xhat_0))*0.1; % initial state estimate covariance matrix
%H = eye(length(xhat_0));      % measurement/observation matrix 
H = [1 0 0;...                 % measurement/observation matrix with constraint
     0 1 0;...
     0 0 1;...
     1 -deltaT/2 deltaT/2];


% Define the constant and error matrices
c_K = [0 ; c_F_in ; c_F_out];  % Constants in AR model
w_K = [0.0000009 ; 296; 319]; % Variance from AR model

% Define the noise
% Measurement noise from instrumentation data sheets
L_sensor_noise    = 0.01;           % std dev of measurement noise, level sensor error
F_in_meter_noise  = 0.01;           % std dev of measurement noise, inlet stream flowmeter error
F_out_meter_noise = 0.01;           % std dev of measurement noise, outlet stream flowmeter error
R = diag([L_sensor_noise^2 ;...
    F_in_meter_noise^2 ;...
    F_out_meter_noise^2 ;...
    0]);                             % measurement noise covariance (measurement error squared)

Q = diag([w_K(1) ; w_K(2); w_K(3)]); % total covariance based on variance from AR model
                                     % process noise covariance 

                                                 
% Define the transition matrix
A = [1 C_L -C_L;...    
     0 a_F_in 0;...
     0 0 a_F_out];


% Define the state estimates vs measurement matrices
x_K = [k_L ; k_F_in ; k_F_out]; % State vector estimate/model prediction

P = P0;  % Setting initial state estimate covariance matrix values 

for n = 2:length(t)-1
   % Priori predictions
   x_K(:,n) = A*x_K(:,n-1) + c_K; % priori predicted state
   P = A.*P.*A' + Q;              % priori predicted covariance

   % Kalman gain
   K = P*H'*inv(H*P*H' + R);

   % Create a set of 'perfect measurements' as a constraint
   Meas_Constraint = x_K(1,n-1) + deltaT/2.*(x_K(2,n-1) - x_K(3,n-1));

   % Redefine the measurement matrix
   z_K(:,n) = [Meas_L(n); Meas_F_in(n) ; Meas_F_out(n) ; Meas_Constraint];

   % Correction
   e(:,n) = z_K(:,n) - (H*x_K(:,n)); % Measurement residual where z is the observation
                                     % vector or measurements
   x_K(:,n) = x_K(:,n) + (K*e(:,n)); % posteriori state estimate
   P = P - (K*H*P);                  % posteriori covariance
end

x_K_L     = x_K(1,:); % Kalman estimate for SD Level
x_K_F_in  = x_K(2,:); % Kalman estimate for SD Inlet Flowrate
x_K_F_out = x_K(3,:); % Kalman estimate for SD Outlet Flowrate

% Plot the results of the Kalman Estimate vs the Raw Data
figure(1)
title('Kalman Filter');

subplot(3,1,1)
plot(t, x_K_L', t, Meas_L, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('L_P_T_A (m)')

subplot(3,1,2)
plot(t, x_K_F_in',t, Meas_F_in, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('F_i_n_P_T (L/s)')

subplot(3,1,3)
plot(t, x_K_F_out', t, Meas_F_out, '--')
legend('Kalman Filter', 'Raw Data')
xlabel('Time (s)')
ylabel('F_o_u_t_P_T (L/s)')

u.L_SD_filtered  = griddedInterpolant(t, x_K_L');
u.F_in_filtered  = griddedInterpolant(t, x_K_F_in');
u.F_out_filtered = griddedInterpolant(t, x_K_F_out');

u.F_in_generated = griddedInterpolant(t, g_F_in');
u.F_out_generated = griddedInterpolant(t, g_F_out');

figure(2)
plot(t, x_K_F_in, t, x_K_F_out);
legend('Fin','Fout'); % Compare the inlet and outlet filtered flowrate values

save KalmanFilterSD.mat u