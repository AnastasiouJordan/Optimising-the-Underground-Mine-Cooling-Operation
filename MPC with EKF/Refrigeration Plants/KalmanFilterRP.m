function x = KalmanFilterRP(x, s, output, u, z, t, p, response)

v = RPIntermediatesKF(x, u, p, t, output, s, response);

% Fridge Plant 1
z_vec1 = [z.T1(end); z.Qrefr1(end); z.Qamb1(end); z.Tin1(end)];        % Redefine the measurement matrix with latest measurements
x01 = [x.T1(end); x.Qrefr1(end); x.Qamb1(end); x.Tin1(end); x.P1(end)]; % Define the starting point for the KF integration
[~,x_integrate1] = ode45(@(t,x) KFFridgePlantsODEs1(s,p,x,u,t,output,1), t, x01);
x_T1 = x_integrate1(end,1);                         % Define the new state estimate for temperature
x_Qrefr1 = x_integrate1(end,2);                     % Define the new state estimate for Q_refr
x_Qamb1 = x_integrate1(end,3);                      % Define the new state estimate for Q_amb
x_Tin1  = v.T_inRP(end,1); 
x_K1 = [x_T1 ; x_Qrefr1; x_Qamb1; x_Tin1];          % State vector estimate/model prediction
P_pri1 = eye(length(p.xhat_0))*x_integrate1(end,5); % Priori predicted covariance
K1 = P_pri1.*p.H'.*inv(p.H.*P_pri1.*p.H' + p.R);    % Kalman gain calculation
e1 = z_vec1 - (p.H.*x_K1);                          % Measurement residual where z is the observation vector or measurements
x_K1 = x_K1 + (K1*e1);                              % Posteriori state estimate
P_post1 = P_pri1 - (K1.*p.H.*P_pri1);               % Posteriori covariance
x.T1(end+1) = x_K1(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr1(end+1) = x_K1(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb1(end+1) = x_K1(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin1(end+1) = x_K1(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P1(end+1) = P_post1(1,1);                         % Kalman estimate assigned to updated state for covariance

% Fridge Plant 2
z_vec2 = [z.T2(end); z.Qrefr2(end); z.Qamb2(end); z.Tin2(end)];        % Redefine the measurement matrix with latest measurements
x02 = [x.T2(end); x.Qrefr2(end); x.Qamb2(end); x.Tin2(end); x.P2(end)]; % Define the starting point for the KF integration
[~,x_integrate2] = ode45(@(t,x) KFFridgePlantsODEs2(s,p,x,u,t,output,2), t, x02);
x_T2 = x_integrate2(end,1);                         % Define the new state estimate for temperature
x_Qrefr2 = x_integrate2(end,2);                     % Define the new state estimate for Q_refr
x_Qamb2 = x_integrate2(end,3);                      % Define the new state estimate for Q_amb
x_Tin2  = v.T_inRP(end,2); 
x_K2 = [x_T2 ; x_Qrefr2; x_Qamb2; x_Tin2];          % State vector estimate/model prediction
P_pri2 = eye(length(p.xhat_0))*x_integrate2(end,3); % Priori predicted covariance
K2 = P_pri2.*p.H'.*inv(p.H.*P_pri2.*p.H' + p.R);    % Kalman gain calculation
e2 = z_vec2 - (p.H.*x_K2);                          % Measurement residual where z is the observation vector or measurements
x_K2 = x_K2 + (K2*e2);                              % Posteriori state estimate
P_post2 = P_pri2 - (K2.*p.H.*P_pri2);               % Posteriori covariance
x.T2(end+1) = x_K2(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr2(end+1) = x_K2(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb2(end+1) = x_K2(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin2(end+1) = x_K2(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P2(end+1) = P_post2(1,1);                         % Kalman estimate assigned to updated state for covariance

% Fridge Plant 3
z_vec3 = [z.T3(end); z.Qrefr3(end); z.Qamb3(end); z.Tin3(end)];        % Redefine the measurement matrix with latest measurements
x03 = [x.T3(end); x.Qrefr3(end); x.Qamb3(end); x.Tin3(end); x.P3(end)]; % Define the starting point for the KF integration
[~,x_integrate3] = ode45(@(t,x) KFFridgePlantsODEs3(s,p,x,u,t,output,3), t, x03);
x_T3 = x_integrate3(end,1);                         % Define the new state estimate for temperature
x_Qrefr3 = x_integrate3(end,2);                     % Define the new state estimate for Q_refr
x_Qamb3 = x_integrate3(end,3);                      % Define the new state estimate for Q_amb
x_Tin3  = v.T_inRP(end,3); 
x_K3 = [x_T3 ; x_Qrefr3; x_Qamb3; x_Tin3];          % State vector estimate/model prediction
P_pri3 = eye(length(p.xhat_0))*x_integrate3(end,3); % Priori predicted covariance
K3 = P_pri3.*p.H'.*inv(p.H.*P_pri3.*p.H' + p.R);    % Kalman gain calculation
e3 = z_vec3 - (p.H.*x_K3);                          % Measurement residual where z is the observation vector or measurements
x_K3 = x_K3 + (K3*e3);                              % Posteriori state estimate
P_post3 = P_pri3 - (K3.*p.H.*P_pri3);               % Posteriori covariance
x.T3(end+1) = x_K3(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr3(end+1) = x_K3(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb3(end+1) = x_K3(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin3(end+1) = x_K3(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P3(end+1) = P_post3(1,1);                         % Kalman estimate assigned to updated state for covariance

% Fridge Plant 4
z_vec4 = [z.T4(end); z.Qrefr4(end); z.Qamb4(end); z.Tin4(end)];        % Redefine the measurement matrix with latest measurements
x04 = [x.T4(end); x.Qrefr4(end); x.Qamb4(end); x.Tin4(end); x.P4(end)]; % Define the starting point for the KF integration
[~,x_integrate4] = ode45(@(t,x) KFFridgePlantsODEs4(s,p,x,u,t,output,4), t, x04);
x_T4 = x_integrate4(end,1);                         % Define the new state estimate for temperature
x_Qrefr4 = x_integrate4(end,2);                     % Define the new state estimate for Q_refr
x_Qamb4 = x_integrate4(end,3);                      % Define the new state estimate for Q_amb
x_Tin4  = v.T_inRP(end,4); 
x_K4 = [x_T4 ; x_Qrefr4; x_Qamb4; x_Tin4];          % State vector estimate/model prediction
P_pri4 = eye(length(p.xhat_0))*x_integrate4(end,3); % Priori predicted covariance
K4 = P_pri4.*p.H'.*inv(p.H.*P_pri4.*p.H' + p.R);    % Kalman gain calculation
e4 = z_vec4 - (p.H.*x_K4);                          % Measurement residual where z is the observation vector or measurements
x_K4 = x_K4 + (K4*e4);                              % Posteriori state estimate
P_post4 = P_pri4 - (K4.*p.H.*P_pri4);               % Posteriori covariance
x.T4(end+1) = x_K4(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr4(end+1) = x_K4(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb4(end+1) = x_K4(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin4(end+1) = x_K4(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P4(end+1) = P_post4(1,1);                         % Kalman estimate assigned to updated state for covariance

% Fridge Plant 5
z_vec5 = [z.T5(end); z.Qrefr5(end); z.Qamb5(end); z.Tin5(end)];        % Redefine the measurement matrix with latest measurements
x05 = [x.T5(end); x.Qrefr5(end); x.Qamb5(end); x.Tin5(end); x.P5(end)]; % Define the starting point for the KF integration
[~,x_integrate5] = ode45(@(t,x) KFFridgePlantsODEs5(s,p,x,u,t,output,5), t, x05);
x_T5 = x_integrate5(end,1);                         % Define the new state estimate for temperature
x_Qrefr5 = x_integrate5(end,2);                     % Define the new state estimate for Q_refr
x_Qamb5 = x_integrate5(end,3);                      % Define the new state estimate for Q_amb
x_Tin5  = v.T_inRP(end,5); 
x_K5 = [x_T5 ; x_Qrefr5; x_Qamb5; x_Tin5];          % State vector estimate/model prediction
P_pri5 = eye(length(p.xhat_0))*x_integrate5(end,3); % Priori predicted covariance
K5 = P_pri5.*p.H'.*inv(p.H.*P_pri5.*p.H' + p.R);    % Kalman gain calculation
e5 = z_vec5 - (p.H.*x_K5);                          % Measurement residual where z is the observation vector or measurements
x_K5 = x_K5 + (K5*e5);                              % Posteriori state estimate
P_post5 = P_pri5 - (K5.*p.H.*P_pri5);               % Posteriori covariance
x.T5(end+1) = x_K5(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr5(end+1) = x_K5(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb5(end+1) = x_K5(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin5(end+1) = x_K5(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P5(end+1) = P_post5(1,1);                         % Kalman estimate assigned to updated state for covariance

end