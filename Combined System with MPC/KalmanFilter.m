function x = KalmanFilter(x, s, output, u, z, t, p, response)

% Square Dam
z_vecSD = [z.L_SD(end)/100*p.height_SD; z.F_in_SD(end)]; % Redefine the measurement matrix with latest measurements
x0SD = [x.L_SD(end); x.F_in_SD(end); x.P_SD(end)];       % Provide starting point for integration with latest states
[~,x_integrateSD] = ode45(@(t,x) KFSquareDamODEs(s, p, x, u, t, output), t, x0SD);
x_L_SD = x_integrateSD(end,1)/100*p.height_SD;                  % Convert previous state estimate to height
x_F_in_SD = x_integrateSD(end,2);
x_K_SD = [x_L_SD ; x_F_in_SD];                                  % State vector estimate/model prediction
x_K_SD = p.A_SD.*x_K_SD + p.c_K_SD;                             % Priori predicted state
P_pri_SD = eye(length(p.xhat_0_SD))*x_integrateSD(end,3);       % Priori predicted covariance
K_SD = P_pri_SD.*p.H_SD'.*inv(p.H_SD.*P_pri_SD.*p.H_SD' + p.R_SD); % Kalman gain
e_SD = z_vecSD - (p.H_SD.*x_K_SD);            % Measurement residual where z is the observation vector or measurements                  
x_K_SD = x_K_SD + (K_SD*e_SD);                % Posteriori state estimate
P_post_SD = P_pri_SD - (K_SD.*p.H_SD.*P_pri_SD); % Posteriori covariance
x.L_SD(end+1) = x_K_SD(1)/p.height_SD*100;    % Kalman estimate for SD Level
x.F_in_SD(end+1) = x_K_SD(2,2);                 % Kalman estimate for SD Inlet Flowrate
x.P_SD(end+1) = P_post_SD(1,1);

% Pre-cooling Towers
z_vecPT = [z.L_PT(end)/100*p.height_PT; z.F_AntiSurge(end); z.T_AntiSurge(end); z.T_outSD(end); z.T_PT(end); z.T_amb(end); z.H_Env(end)]; % Redefine the measurement matrix with latest measurements
x0PT = [x.L_PT(end); x.F_AntiSurge(end); x.T_AntiSurge(end); x.T_outSD(end); x.T_PT(end); x.T_amb(end); x.H_Env(end); x.P_PT(end)];      % Provide starting point for integration with latest states
[~,x_integratePT] = ode45(@(t,x) KFPreCoolingTowerODEs(s, p, x, u, t, output), t, x0PT);
x_L_PT = x_integratePT(end,1)/100*p.height_PT; % Convert previous state estimate to height
x_F_AntiSurge = x_integratePT(end,2);
x_T_AntiSurge = x_integratePT(end,3);
x_T_outSD = x_integratePT(end,4);
x_T_PT    = x_integratePT(end,5);
x_T_amb   = x_integratePT(end,6);
x_H_Env   = x_integratePT(end,7);
x_K_PT = [x_L_PT ; x_F_AntiSurge; x_T_AntiSurge; x_T_outSD; x_T_PT; x_T_amb; x_H_Env]; % State vector estimate/model prediction
P_pri_PT = eye(length(p.xhat_0_PT))*x_integratePT(end,8);            % Priori predicted covariance
K_PT = P_pri_PT.*p.H_PT'.*inv(p.H_PT.*P_pri_PT.*p.H_PT' + p.R_PT);   % Kalman gain
e_PT = z_vecPT - (p.H_PT.*x_K_PT);               % Measurement residual where z is the observation vector or measurements
x_K_PT = x_K_PT + (K_PT*e_PT);                   % Posteriori state estimate
P_post_PT = P_pri_PT - (K_PT.*p.H_PT.*P_pri_PT); % Posteriori covariance
x.L_PT(end+1) = x_K_PT(1,1)/p.height_PT*100;       % Kalman estimate for PT Level
x.F_AntiSurge(end+1) = x_K_PT(2,2);  
x.T_AntiSurge(end+1) = x_K_PT(3,3);
x.T_outSD(end+1) = x_K_PT(4,4);
x.T_PT(end+1) = x_K_PT(5,5);
x.T_amb(end+1) = x_K_PT(6,6);
x.H_Env(end+1) = x_K_PT(7,7);
x.P_PT(end+1) = P_post_PT(1,1);



% Fridge Plants
v = RPIntermediatesKF(x, u, p, t, output, s, response);
% Fridge Plant 1
z_vec1 = [z.T1(end); z.Qrefr1(end); z.Qamb1(end); z.Tin1(end)];         % Redefine the measurement matrix with latest measurements
x01 = [x.T1(end); x.Qrefr1(end); x.Qamb1(end); x.Tin1(end); x.P1(end)]; % Define the starting point for the KF integration
[~,x_integrate1] = ode45(@(t,x) KFFridgePlantsODEs1(s,p,x,u,t,output,1,v), t, x01);
x_T1 = x_integrate1(end,1);                         % Define the new state estimate for temperature
x_Qrefr1 = x_integrate1(end,2);                     % Define the new state estimate for Q_refr
x_Qamb1 = x_integrate1(end,3);                      % Define the new state estimate for Q_amb
x_Tin1  = v.T_inRP(end,1); 
x_K1 = [x_T1 ; x_Qrefr1; x_Qamb1; x_Tin1];          % State vector estimate/model prediction
P_pri1 = eye(length(p.xhat_0_RP))*x_integrate1(end,5); % Priori predicted covariance
K1 = P_pri1.*p.H_RP'.*inv(p.H_RP.*P_pri1.*p.H_RP' + p.R_RP);    % Kalman gain calculation
e1 = z_vec1 - (p.H_RP.*x_K1);                          % Measurement residual where z is the observation vector or measurements
x_K1 = x_K1 + (K1*e1);                              % Posteriori state estimate
P_post1 = P_pri1 - (K1.*p.H_RP.*P_pri1);               % Posteriori covariance
x.T1(end+1) = x_K1(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr1(end+1) = x_K1(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb1(end+1) = x_K1(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin1(end+1) = x_K1(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P1(end+1) = P_post1(1,1);                         % Kalman estimate assigned to updated state for covariance

% Fridge Plant 2
z_vec2 = [z.T2(end); z.Qrefr2(end); z.Qamb2(end); z.Tin2(end)];        % Redefine the measurement matrix with latest measurements
x02 = [x.T2(end); x.Qrefr2(end); x.Qamb2(end); x.Tin2(end); x.P2(end)]; % Define the starting point for the KF integration
[~,x_integrate2] = ode45(@(t,x) KFFridgePlantsODEs2(s,p,x,u,t,output,2,v), t, x02);
x_T2 = x_integrate2(end,1);                         % Define the new state estimate for temperature
x_Qrefr2 = x_integrate2(end,2);                     % Define the new state estimate for Q_refr
x_Qamb2 = x_integrate2(end,3);                      % Define the new state estimate for Q_amb
x_Tin2  = v.T_inRP(end,2); 
x_K2 = [x_T2 ; x_Qrefr2; x_Qamb2; x_Tin2];          % State vector estimate/model prediction
P_pri2 = eye(length(p.xhat_0_RP))*x_integrate2(end,3); % Priori predicted covariance
K2 = P_pri2.*p.H_RP'.*inv(p.H_RP.*P_pri2.*p.H_RP' + p.R_RP);    % Kalman gain calculation
e2 = z_vec2 - (p.H_RP.*x_K2);                          % Measurement residual where z is the observation vector or measurements
x_K2 = x_K2 + (K2*e2);                              % Posteriori state estimate
P_post2 = P_pri2 - (K2.*p.H_RP.*P_pri2);               % Posteriori covariance
x.T2(end+1) = x_K2(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr2(end+1) = x_K2(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb2(end+1) = x_K2(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin2(end+1) = x_K2(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P2(end+1) = P_post2(1,1);                         % Kalman estimate assigned to updated state for covariance

% Fridge Plant 3
z_vec3 = [z.T3(end); z.Qrefr3(end); z.Qamb3(end); z.Tin3(end)];         % Redefine the measurement matrix with latest measurements
x03 = [x.T3(end); x.Qrefr3(end); x.Qamb3(end); x.Tin3(end); x.P3(end)]; % Define the starting point for the KF integration
[~,x_integrate3] = ode45(@(t,x) KFFridgePlantsODEs3(s,p,x,u,t,output,3,v), t, x03);
x_T3 = x_integrate3(end,1);                         % Define the new state estimate for temperature
x_Qrefr3 = x_integrate3(end,2);                     % Define the new state estimate for Q_refr
x_Qamb3 = x_integrate3(end,3);                      % Define the new state estimate for Q_amb
x_Tin3  = v.T_inRP(end,3); 
x_K3 = [x_T3 ; x_Qrefr3; x_Qamb3; x_Tin3];          % State vector estimate/model prediction
P_pri3 = eye(length(p.xhat_0_RP))*x_integrate3(end,3); % Priori predicted covariance
K3 = P_pri3.*p.H_RP'.*inv(p.H_RP.*P_pri3.*p.H_RP' + p.R_RP);    % Kalman gain calculation
e3 = z_vec3 - (p.H_RP.*x_K3);                          % Measurement residual where z is the observation vector or measurements
x_K3 = x_K3 + (K3*e3);                              % Posteriori state estimate
P_post3 = P_pri3 - (K3.*p.H_RP.*P_pri3);               % Posteriori covariance
x.T3(end+1) = x_K3(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr3(end+1) = x_K3(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb3(end+1) = x_K3(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin3(end+1) = x_K3(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P3(end+1) = P_post3(1,1);                         % Kalman estimate assigned to updated state for covariance

% Fridge Plant 4
z_vec4 = [z.T4(end); z.Qrefr4(end); z.Qamb4(end); z.Tin4(end)];         % Redefine the measurement matrix with latest measurements
x04 = [x.T4(end); x.Qrefr4(end); x.Qamb4(end); x.Tin4(end); x.P4(end)]; % Define the starting point for the KF integration
[~,x_integrate4] = ode45(@(t,x) KFFridgePlantsODEs4(s,p,x,u,t,output,4,v), t, x04);
x_T4 = x_integrate4(end,1);                         % Define the new state estimate for temperature
x_Qrefr4 = x_integrate4(end,2);                     % Define the new state estimate for Q_refr
x_Qamb4 = x_integrate4(end,3);                      % Define the new state estimate for Q_amb
x_Tin4  = v.T_inRP(end,4); 
x_K4 = [x_T4 ; x_Qrefr4; x_Qamb4; x_Tin4];          % State vector estimate/model prediction
P_pri4 = eye(length(p.xhat_0_RP))*x_integrate4(end,3); % Priori predicted covariance
K4 = P_pri4.*p.H_RP'.*inv(p.H_RP.*P_pri4.*p.H_RP' + p.R_RP);    % Kalman gain calculation
e4 = z_vec4 - (p.H_RP.*x_K4);                          % Measurement residual where z is the observation vector or measurements
x_K4 = x_K4 + (K4*e4);                              % Posteriori state estimate
P_post4 = P_pri4 - (K4.*p.H_RP.*P_pri4);               % Posteriori covariance
x.T4(end+1) = x_K4(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr4(end+1) = x_K4(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb4(end+1) = x_K4(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin4(end+1) = x_K4(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P4(end+1) = P_post4(1,1);                         % Kalman estimate assigned to updated state for covariance

% Fridge Plant 5
z_vec5 = [z.T5(end); z.Qrefr5(end); z.Qamb5(end); z.Tin5(end)];         % Redefine the measurement matrix with latest measurements
x05 = [x.T5(end); x.Qrefr5(end); x.Qamb5(end); x.Tin5(end); x.P5(end)]; % Define the starting point for the KF integration
[~,x_integrate5] = ode45(@(t,x) KFFridgePlantsODEs5(s,p,x,u,t,output,5,v), t, x05);
x_T5 = x_integrate5(end,1);                         % Define the new state estimate for temperature
x_Qrefr5 = x_integrate5(end,2);                     % Define the new state estimate for Q_refr
x_Qamb5 = x_integrate5(end,3);                      % Define the new state estimate for Q_amb
x_Tin5  = v.T_inRP(end,5); 
x_K5 = [x_T5 ; x_Qrefr5; x_Qamb5; x_Tin5];          % State vector estimate/model prediction
P_pri5 = eye(length(p.xhat_0_RP))*x_integrate5(end,3); % Priori predicted covariance
K5 = P_pri5.*p.H_RP'.*inv(p.H_RP.*P_pri5.*p.H_RP' + p.R_RP);    % Kalman gain calculation
e5 = z_vec5 - (p.H_RP.*x_K5);                          % Measurement residual where z is the observation vector or measurements
x_K5 = x_K5 + (K5*e5);                              % Posteriori state estimate
P_post5 = P_pri5 - (K5.*p.H_RP.*P_pri5);               % Posteriori covariance
x.T5(end+1) = x_K5(1);                              % Kalman estimate assigned to updated state for temperature
x.Qrefr5(end+1) = x_K5(2,2);                        % Kalman estimate assigned to updated state for Q_refr
x.Qamb5(end+1) = x_K5(3,3);                         % Kalman estimate assigned to updated state for Q_amb
x.Tin5(end+1) = x_K5(4,4);                          % Kalman estimate assigned to updated state for T_in
x.P5(end+1) = P_post5(1,1);                         % Kalman estimate assigned to updated state for covariance

% Chill Dam
v.T_PT = x.T_PT(end); v.T1 = x.T1(end); v.T2 = x.T2(end); v.T3 = x.T3(end); 
v.T4 = x.T4(end); v.T5 = x.T5(end); v.T_amb = x.T_amb(end);
z_vecCD = [z.L_CD(end)/100*p.height_CD; z.T_CD(end); z.F_Ice(end); z.F_UG(end)]; % Redefine the measurement matrix with latest measurements
x0CD = [x.L_CD(end); x.T_CD(end); x.F_Ice(end); x.F_UG(end); x.P_CD(end)]; % Provide starting point for integration with latest states
[~,x_integrateCD] = ode45(@(t,x) KFChillDamODEs(s, p, x, u, t, output, response, v), t, x0CD);
x_L_CD = x_integrateCD(end,1)/100*p.height_CD;                  % Convert previous state estimate to height
x_T_CD = x_integrateCD(end,2);
x_F_Ice = x_integrateCD(end,3);
x_F_UG = x_integrateCD(end,4);
x_K_CD = [x_L_CD ; x_T_CD; x_F_Ice; x_F_UG];                            % State vector estimate/model prediction
P_pri_CD = eye(length(p.xhat_0_CD))*x_integrateCD(end,5);       % Priori predicted covariance
K_CD = P_pri_CD.*p.H_CD'.*inv(p.H_CD.*P_pri_CD.*p.H_CD' + p.R_CD);
e_CD = z_vecCD - (p.H_CD.*x_K_CD);                              % Measurement residual where z is the observation vector or measurements
x_K_CD = x_K_CD + (K_CD*e_CD);                                  % Posteriori state estimate
P_post_CD = P_pri_CD - (K_CD.*p.H_CD.*P_pri_CD);                % Posteriori covariance
x.L_CD(end+1) = x_K_CD(1,1)/p.height_CD*100;                      % Kalman estimate for CD Level
x.T_CD(end+1) = x_K_CD(2,2);
x.F_Ice(end+1) = x_K_CD(3,3);                                   % Kalman estimate for CD Inlet Flowrate
x.F_UG(end+1) = x_K_CD(4,4);
x.P_CD(end+1) = P_post_CD(1,1);

end