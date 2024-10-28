function x = KalmanFilterPT(x, s, output, u, z, t, p)

% Define the inital conditions
% States: x = [x1; x2; x3]
% x1 = L
% x2 = Fin
% x3 = Fout

z_vec = [z.L(end)/100*p.height_CD; z.F_in(end)]; % Redefine the measurement matrix with latest measurements

x0 = [x.L(end); x.F_in(end); x.h_CD(end); x.P(end)];
%x0 = [x.L(end); x.F_in(end); x.P(end)];

[~,x_integrate] = ode45(@(t,x) KFChillDamODEs(s, p, x, u, t, output), t, x0);


% Define the state estimates vs measurement matrices
x_L = x_integrate(end,1)/100*p.height_CD; % Convert previous state estimate to height
x_F_in = x_integrate(end,2);
x_K = [x_L ; x_F_in]; % State vector estimate/model prediction

% For loop no longer needed, but it would have started here
   % Priori predictions
   x_K = p.A.*x_K + p.c_K; % priori predicted state
   %P_pri = p.A.*eye(length(p.xhat_0))*x.P(end).*p.A' + p.Q;
   P_pri = eye(length(p.xhat_0))*x_integrate(end,4);            % priori predicted covariance

   % Kalman gain
   K = P_pri.*p.H'.*inv(p.H.*P_pri.*p.H' + p.R);

   % Create a set of 'perfect measurements' as a constraint
   %Meas_Constraint = x.K(1,n-1) + deltaT/2.*(x.K(2,n-1) - x.K(3,n-1));
   

   % Correction
   e = z_vec - (p.H.*x_K); % Measurement residual where z is the observation
                                     % vector or measurements
   x_K = x_K + (K*e); % posteriori state estimate
   P_post = P_pri - (K.*p.H.*P_pri);                  % posteriori covariance

x.L(end+1) = x_K(1)/p.height_CD*100; % Kalman estimate for PT Level
x.F_in(end+1) = x_K(2); % Kalman estimate for SD Inlet Flowrate
x.P(end+1) = P_post(1,1);
%x.h_PT(end+1) = x_integrate(end,3);


% % Plot the results of the Kalman Estimate vs the Raw Data
% figure(1)
% title('Kalman Filter');
% 
% subplot(3,1,1)
% plot(t, x.K_L', t, Meas_L, '--')
% legend('Kalman Filter', 'Raw Data')
% xlabel('Time (s)')
% ylabel('L_P_T_A (m)')
% 
% subplot(3,1,2)
% plot(t, x.K_F_in',t, Meas_F_in, '--')
% legend('Kalman Filter', 'Raw Data')
% xlabel('Time (s)')
% ylabel('F_i_n_P_T (L/s)')
% 
% subplot(3,1,3)
% plot(t, x.K_F_out', t, Meas_F_out, '--')
% legend('Kalman Filter', 'Raw Data')
% xlabel('Time (s)')
% ylabel('F_o_u_t_P_T (L/s)')
% 
% u.L_SD_filtered  = griddedInterpolant(t, x.K_L');
% u.F_in_filtered  = griddedInterpolant(t, x.K_F_in');
% u.F_out_filtered = griddedInterpolant(t, x.K_F_out');
% 
% u.F_in_generated = griddedInterpolant(t, g.F_in');
% u.F_out_generated = griddedInterpolant(t, g.F_out');
% 
% figure(2)
% plot(t, x.K_F_in, t, x.K_F_out);
% legend('Fin','Fout'); % Compare the inlet and outlet filtered flowrate values
% 
% save KalmanFilterSD.mat u
end