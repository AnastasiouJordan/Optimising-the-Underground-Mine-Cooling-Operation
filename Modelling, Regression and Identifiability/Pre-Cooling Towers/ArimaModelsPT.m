clc
clear

load SavedInterpolantsPT.mat

% ARIMA Models

% Define Dimensions
p.height_PT = 3;     % m, Height of the Pre-cooling Tower basins
p.length_PT = 14.8;  % m, Length of the Pre-cooling Tower basins
p.width_PT  = 10.45*2;  % m, Width of the Pre-cooling Tower basins
p.area_PT   = p.length_PT*p.width_PT; % m2, Area of the Pre-cooling Tower basins
p.volume_PT = p.area_PT*p.height_PT;  % m3, Volume of the Pre-cooling Tower basins


% ARIMA MODEL FOR F_inPT
% Load the data
Raw_F_in  = u.F_outSD(t) + u.F_AntiSurge(t);
Size_F_in = size(Raw_F_in,1);
t1 = t(1:51086,1);

% x = Raw_F_in;
% windowSize = 1000; 
% b = (1/windowSize)*ones(1,windowSize);
% a = 1;
% y = filter(b,a,x);
% [d1,p1] = aryule(y,1);
% 
% plot(t, y,'LineWidth', 2);
% hold on
% plot(t, x);
% hold off


% Create Model Template
Mdl_F_in = arima(1,0,0); % First-order Auto-Regressive Model plus constant

% Here we specify the following arguments for the model:
% p-value: AR - the autoregressive polynomial degreee
% D-value: MA - the moving average degree of integration
% q-value: ARMA - the moving average polynomial degree

% Here, we have:
% 1 nonseasonal AR polynomial lag
% 0 degree nonseasonal integration polynomial
% 0 nonseasonal MA polynomial lags


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_F_in = 1:Mdl_F_in.P;               % Contains 2 observations
estsample_F_in = (Mdl_F_in.P + 1):Size_F_in; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_F_in = estimate(Mdl_F_in, Raw_F_in(estsample_F_in), 'Y0', Raw_F_in(presample_F_in));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
a_F_in = cell2mat(EstMdl_F_in.AR(1)); % AR coefficient, alpha
w_F_in = EstMdl_F_in.Variance;        % Variance which is rooted to obtain error
c_F_in = EstMdl_F_in.Constant;        % Constant or intercept


% ARIMA MODEL FOR F_outPT

% Load the data
Raw_F_out = u.F_outPT(t);
Size_F_out = size(Raw_F_out,1);

% Create Model Template
Mdl_F_out = arima(1,0,0); % First-order Auto-Regressive Model plus constant

% Here we specify the following arguments for the model:
% p-value: AR - the autoregressive polynomial degreee
% D-value: MA - the moving average degree of integration
% q-value: ARMA - the moving average polynomial degree

% Here, we have:
% 1 nonseasonal AR polynomial lag
% 0 degree nonseasonal integration polynomial
% 0 nonseasonal MA polynomial lags


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_F_out = 1:Mdl_F_out.P;                % Contains 2 observations
estsample_F_out = (Mdl_F_out.P + 1):Size_F_out; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_F_out = estimate(Mdl_F_out, Raw_F_out(estsample_F_out), 'Y0', Raw_F_out(presample_F_out));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
a_F_out = cell2mat(EstMdl_F_out.AR(1)); % AR coefficient, alpha
w_F_out = EstMdl_F_out.Variance;        % Variance which is rooted to obtain error
c_F_out = EstMdl_F_out.Constant;        % Constant or intercept


% ARIMA MODEL FOR L_PT (just using an average of the two levels)

% Load the data
Raw_L = (u.L_PTA(t) + u.L_PTB(t))./2; % Raw level measurements as percentages
Raw_L = Raw_L/100*p.height_PT;        % Raw level measurements as heighs in meters
Size_L = size(Raw_L,1);

% Create Model Template
dep = 0.8; % Dependence of the next L measurement on the previous one
%Mdl_L_A = arima(1,0,0); % This gives an alpha value of 0.9999 and poor fit
Mdl_L = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant

% Here we specify the following arguments for the model:
% p-value: AR - the autoregressive polynomial degreee
% D-value: MA - the moving average degree of integration
% q-value: ARMA - the moving average polynomial degree

% Here, we have:
% 1 nonseasonal AR polynomial lag
% 0 degree nonseasonal integration polynomial
% 0 nonseasonal MA polynomial lags
% We also specify the alpha value as 0.99


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_L = 1:Mdl_L.P;            % Contains 2 observations
estsample_L = (Mdl_L.P + 1):Size_L; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_L = estimate(Mdl_L, Raw_L(estsample_L), 'Y0', Raw_L(presample_L));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
a_L    = cell2mat(EstMdl_L.AR(1)); % AR coefficient, alpha
w_L    = EstMdl_L.Variance;        % Variance which is rooted to obtain error
c_L    = EstMdl_L.Constant;        % Constant or intercept
deltaT = t(2) - t(1);                % s, Delta time in seconds
C_L    = (deltaT/p.area_PT)*0.001;   % m, Coefficient to the flowrates.
                                     % This calculation takes delta time and
                                     % divides it by the area of the Chill Dam
                                     % in m2. This is then multiplied by 0.001
                                     % to adjust for the flowrate which is in
                                     % L/s. This gives a constant in meters,
                                     % equivalent to the units in which level
                                     % is measured.


% AR Model Predictions
k_F_in(1)  = Raw_F_in(1);  % Specify the inital value of the inlet flowrate
k_F_out(1) = Raw_F_out(1); % Specify the inital value of the outlet flowrate
k_L(1)     = Raw_L(1);     % Specify the inital value of the level
for i = 1:51086            % From here onwards, make predictions
    e_F_in(i)    = randn*sqrt(w_F_in);
    k_F_in(i+1)  = c_F_in + a_F_in*k_F_in(i) + e_F_in(i);
    e_F_out(i)   = randn*sqrt(w_F_out);
    k_F_out(i+1) = c_F_out + a_F_out*k_F_out(i) + e_F_out(i);
    e_L(i)       = randn*sqrt(w_L);
    k_L(i+1)     = a_L*k_L(i) + C_L*(k_F_in(i) - k_F_out(i)) + e_L(i) + c_L;
end
g_F_in(1)  = Raw_F_in(1);  % Specify the inital value of the inlet flowrate
g_F_out(1) = Raw_F_out(1); % Specify the inital value of the outlet flowrate
g_L(1)     = Raw_L(1);     % Specify the inital value of the level
for i = 1:51086           
    e_F_in(i)    = randn*sqrt(w_F_in);
    g_F_in(i+1)  = c_F_in + a_F_in*k_F_in(i) + e_F_in(i);    % Generated data used to test the model in simulation
    e_F_out(i)   = randn*sqrt(w_F_out);
    g_F_out(i+1) = c_F_out + a_F_out*k_F_out(i) + e_F_out(i); % Generated data used to test the model in simulation
    e_L(i)       = randn*sqrt(w_L);
    g_L(i+1)     = a_L*k_L(i) + C_L*(k_F_in(i) - k_F_out(i)) + e_L(i) + c_L; % Generated data used to test the model in simulation
end
% Where in the above loop:
% e(i) is the error obtained from the variance calculated by the model for each variable.
% c is the constant or intercept calculated by the model for each variable.
% a is the 'alpha' co-efficient for the dependence of the next value on the previous one for each variable.
% k(i) is the prediction for each variable.
% C is the constant used only in the level predictions, derived from the relationship between level and flowrates.



% INLET FLOWRATE
% Plot Actual Data vs AR Model
figure(1)
title('Inlet Flowrate AR Model');
subplot(3,1,1)
plot(t, Raw_F_in, 'r', t, k_F_in', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('F_i_n_P_T (L/s)')

% Plot Observations and Fitted Values
subplot(3,1,2)
resid_F_in = infer(EstMdl_F_in, Raw_F_in(estsample_F_in), ...
             'Y0', Raw_F_in(presample_F_in));      % Infer residuals from
                                                   % the estimated model
yhat_F_in = Raw_F_in(estsample_F_in) - resid_F_in; % Compute the fitted values
plot(t1, Raw_F_in(estsample_F_in),'r', t1, yhat_F_in,'b--','LineWidth',1)
legend('Observations', 'Fitted Values')

% Plot Residuals vs Fitted Values
subplot(3,1,3)
plot(yhat_F_in,resid_F_in,'b.')
ylabel('Residuals')
xlabel('Fitted Values')

% OUTLET FLOWRATE
% Plot Actual Data vs AR Model
figure(2)
title('Outlet Flowrate AR Model');
subplot(3,1,1)
plot(t, Raw_F_out, 'r', t, k_F_out', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('F_o_u_t_P_T (L/s)')

% Plot Observations and Fitted Values
subplot(3,1,2)
resid_F_out = infer(EstMdl_F_out, Raw_F_out(estsample_F_out), ...
             'Y0', Raw_F_out(presample_F_out));        % Infer residuals from
                                                       % the estimated model
yhat_F_out = Raw_F_out(estsample_F_out) - resid_F_out; % Compute the fitted values
plot(t1, Raw_F_out(estsample_F_out),'r', t1, yhat_F_out,'b--','LineWidth',1)
legend('Observations', 'Fitted Values')

% Plot Residuals vs Fitted Values
subplot(3,1,3)
plot(yhat_F_out,resid_F_out,'b.')
ylabel('Residuals')
xlabel('Fitted Values')

% LEVEL
% Plot Actual Data vs AR Model
figure(3)
title('Level AR Model');
subplot(3,1,1)
plot(t, Raw_L, 'r', t, k_L', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('L_P_T (m)')

% Plot Observations and Fitted Values
subplot(3,1,2)
resid_L = infer(EstMdl_L, Raw_L(estsample_L), ...
             'Y0', Raw_L(presample_L)); % Infer residuals from
                                        % the estimated model
yhat_L = Raw_L(estsample_L) - resid_L;  % Compute the fitted values
plot(t1, Raw_L(estsample_L),'r', t1, yhat_L,'b--','LineWidth',1)
legend('Observations', 'Fitted Values')

% Plot Residuals vs Fitted Values
subplot(3,1,3)
plot(yhat_L, resid_L,'b.')
ylabel('Residuals')
xlabel('Fitted Values')


% Convert Level back to %
figure(5)
title('Level AR Model in %');
k_L_percent     = k_L*100/p.height_PT;
L_Data_percent  = Raw_L*100/p.height_PT;
plot(t, L_Data_percent, 'r', t, k_L_percent', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('L_P_T (%)')


save ArimaModelsPT.mat c_F_in c_F_out c_L w_F_in w_F_out w_L C_L a_F_in a_F_out a_L k_F_in k_F_out k_L Raw_F_in Raw_F_out Raw_L deltaT g_F_in g_F_out g_L
