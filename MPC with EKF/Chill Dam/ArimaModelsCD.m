clc
clear

load SavedInterpolantsCD.mat

% ARIMA Models

% Define Dimensions
p.height_CD = 6.7; % m, Height of the Chill Dam
p.length_CD = 48;  % m, Length of the Chill Dam
p.width_CD  = 25;  % m, Width of the Chill Dam
p.area_CD   = p.length_CD*p.width_CD; % m2, Area of the Chill Dam

% ARIMA MODEL FOR F_inCD

% Load the data
Raw_F_in  = u.F_inRP(t) - u.F_UG(t);
Size_F_in = size(Raw_F_in,1);
t1 = t(1:51086,1);

% Create Model Template
Mdl_F_in = arima(1,0,0); % First-order Auto-Regressive Model plus constant

% Here we specify the following arguments for the model:
% p-value: AR - the autoregressive polynomial degree
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


% ARIMA MODEL FOR F_outCD

% Load the data
Raw_F_out = u.F_outCD;
Raw_F_out = Raw_F_out.Values;
Size_F_out = size(Raw_F_out,1);

% Create Model Template
Mdl_F_out = arima(1,0,0); % First-order Auto-Regressive Model plus constant


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


% ARIMA MODEL FOR F_outPT

% Load the data
Raw_F_outPT = u.F_outPT(t);
Size_F_outPT = size(Raw_F_outPT,1);

% Create Model Template
Mdl_F_outPT = arima(1,0,0); % First-order Auto-Regressive Model plus constant


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_F_outPT = 1:Mdl_F_outPT.P;                % Contains 2 observations
estsample_F_outPT = (Mdl_F_outPT.P + 1):Size_F_outPT; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_F_outPT = estimate(Mdl_F_outPT, Raw_F_outPT(estsample_F_outPT), 'Y0', Raw_F_outPT(presample_F_outPT));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
a_F_outPT = cell2mat(EstMdl_F_outPT.AR(1)); % AR coefficient, alpha
w_F_outPT = EstMdl_F_outPT.Variance;        % Variance which is rooted to obtain error
c_F_outPT = EstMdl_F_outPT.Constant;        % Constant or intercept



% ARIMA MODEL FOR F_UG

% Load the data
Raw_F_UG = u.F_UG(t);
Size_F_UG = size(Raw_F_UG,1);

% Create Model Template
Mdl_F_UG = arima(1,0,0); % First-order Auto-Regressive Model plus constant


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_F_UG = 1:Mdl_F_UG.P;                % Contains 2 observations
estsample_F_UG = (Mdl_F_UG.P + 1):Size_F_UG; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_F_UG = estimate(Mdl_F_UG, Raw_F_UG(estsample_F_UG), 'Y0', Raw_F_UG(presample_F_UG));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
a_F_UG = cell2mat(EstMdl_F_UG.AR(1)); % AR coefficient, alpha
w_F_UG = EstMdl_F_UG.Variance;        % Variance which is rooted to obtain error
c_F_UG = EstMdl_F_UG.Constant;        % Constant or intercept



% ARIMA MODEL FOR L_CD

% Load the data
Raw_L = u.L_CD;
Raw_L = Raw_L.Values;          % Raw level measurements as percentages
Raw_L = Raw_L/100*p.height_CD; % Raw level measurements as heighs in meters
Size_L = size(Raw_L,1);

% Create Model Template
dep_L = 0.99; % Dependence of the next L measurement on the previous one
%Mdl_L = arima(1,0,0); % This gives an alpha value of 0.9999 and poor fit
Mdl_L = arima('AR',{dep_L}); % First-order Auto-Regressive Model plus constant


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
a_L = cell2mat(EstMdl_L.AR(1)); % AR coefficient, alpha
w_L = EstMdl_L.Variance;        % Variance which is rooted to obtain error
c_L = EstMdl_L.Constant;        % Constant or intercept
deltaT = t(2) - t(1);           % s, Delta time in seconds
C_L = deltaT/p.area_CD*0.001;   % m, Coefficient to the flowrates.
                                % This calculation takes delta time and
                                % divides it by the area of the Chill Dam
                                % in m2. This is then multiplied by 0.001
                                % to adjust for the flowrate which is in
                                % L/s. This gives a constant in meters,
                                % equivalent to the units in which level
                                % is measured.

% ARIMA MODEL FOR F_Ice

% Load the data

delta_L = diff(Raw_L); 
delta_L(end+1,:) = 0;
Raw_F_Ice = u.F_outPT(t) - u.F_UG(t) - (delta_L.*p.area_CD.*1000./(t(2)-t(1)));
Size_F_Ice = size(Raw_F_Ice,1);
ipt = findchangepts(Raw_F_Ice,'Statistic','std','MinThreshold',10);
intvls = [1 ipt' length(Raw_F_Ice)];
for k1 = 1:length(intvls)-2
    interval_means(k1) = mean(Raw_F_Ice(intvls(k1):intvls(k1+1)-1));
end
interval_means(k1+1) = mean(Raw_F_Ice(intvls(k1):intvls(end)));
xval = [];
yval = [];
for k1 = 1:length(intvls)-1
    xval = intvls(k1):intvls(k1+1);
ynext = interval_means(k1)*ones(1,(length(xval)-1));
yval = [yval ynext];
end
yval(end+1) = 0;

figure(1)
plot(t,Raw_F_Ice);
hold on
plot(t,yval','LineWidth',2);
hold off

Raw_F_Ice = yval';

% Create Model Template
%Mdl_F_Ice = arima(1,0,0); % First-order Auto-Regressive Model plus constant
dep_F_Ice = 0.98; % Dependence of the next L measurement on the previous one
Mdl_F_Ice = arima('AR',{dep_F_Ice}); % First-order Auto-Regressive Model plus constant

% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_F_Ice = 1:Mdl_F_Ice.P;                % Contains 2 observations
estsample_F_Ice = (Mdl_F_Ice.P + 1):Size_F_Ice; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_F_Ice = estimate(Mdl_F_Ice, Raw_F_Ice(estsample_F_Ice), 'Y0', Raw_F_Ice(presample_F_Ice));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
a_F_Ice = cell2mat(EstMdl_F_Ice.AR(1)); % AR coefficient, alpha
w_F_Ice = EstMdl_F_Ice.Variance;        % Variance which is rooted to obtain error
c_F_Ice = EstMdl_F_Ice.Constant;        % Constant or intercept

% AR Model Predictions
% OVER THE CHILL DAM
k_F_in(1)  = Raw_F_in(1);  % Specify the inital value of the inlet flowrate
k_F_out(1) = Raw_F_out(1); % Specify the inital value of the outlet flowrate
k_L1(1)    = Raw_L(1);     % Specify the inital value of the level
for i = 1:51086            % From here onwards, make predictions
    e_F_in(i)    = randn*sqrt(w_F_in);
    k_F_in(i+1)  = c_F_in + a_F_in*k_F_in(i) + e_F_in(i);
    e_F_out(i)   = randn*sqrt(w_F_out);
    k_F_out(i+1) = c_F_out + a_F_out*k_F_out(i) + e_F_out(i);
    e_L(i)       = randn*sqrt(w_L);
    k_L1(i+1)    = a_L*k_L1(i) + C_L*(k_F_in(i) - k_F_out(i)) + e_L(i) + c_L;
end

g_F_in(1)  = Raw_F_in(1);  % Specify the inital value of the inlet flowrate
g_F_out(1) = Raw_F_out(1); % Specify the inital value of the outlet flowrate
g_L1(1)    = Raw_L(1);     % Specify the inital value of the level
for i = 1:51086            % From here onwards, make predictions
    e_F_in(i)    = randn*sqrt(w_F_in);
    g_F_in(i+1)  = c_F_in + a_F_in*k_F_in(i) + e_F_in(i);
    e_F_out(i)   = randn*sqrt(w_F_out);
    g_F_out(i+1) = c_F_out + a_F_out*k_F_out(i) + e_F_out(i);
    e_L(i)       = randn*sqrt(w_L);
    g_L1(i+1)    = a_L*k_L1(i) + C_L*(k_F_in(i) - k_F_out(i)) + e_L(i) + c_L;
end


% Where in the above loop:
% e(i) is the error obtained from the variance calculated by the model for each variable.
% c is the constant or intercept calculated by the model for each variable.
% a is the 'alpha' co-efficient for the dependence of the next value on the previous one for each variable.
% k(i) is the prediction for each variable.
% C is the constant used only in the level predictions, derived from the relationship between level and flowrates.
x_K1 = [k_L1 ; k_F_in ; k_F_out]; % State vector estimate/model prediction

% OVER PT OUTLET, ICE PLANT FLOW AND UG FLOW
k_F_outPT(1) = Raw_F_outPT(1);  % Specify the inital value of the inlet flowrate
k_F_Ice(1)   = Raw_F_Ice(1); % Specify the inital value of the outlet flowrate
k_F_UG(1)    = Raw_F_UG(1);
k_L2(1)      = Raw_L(1);     % Specify the inital value of the level
for i = 1:51086            % From here onwards, make predictions
    e_F_outPT(i)   = randn*sqrt(w_F_outPT);
    k_F_outPT(i+1) = c_F_outPT + a_F_outPT*k_F_outPT(i) + e_F_outPT(i);
    e_F_Ice(i)     = randn*sqrt(w_F_Ice/100);
    k_F_Ice(i+1)   = c_F_Ice + a_F_Ice*k_F_Ice(i) + e_F_Ice(i);
    e_F_UG(i)      = randn*sqrt(w_F_UG);
    k_F_UG(i+1)    = c_F_UG + a_F_UG*k_F_UG(i) + e_F_UG(i);
    e_L(i)         = randn*sqrt(w_L);
    k_L2(i+1)      = a_L*k_L2(i) + C_L*(k_F_outPT(i) - k_F_Ice(i) - k_F_UG(i)) + e_L(i) + c_L;
end

g_F_outPT(1) = Raw_F_outPT(1);  % Specify the inital value of the inlet flowrate
g_F_Ice(1)   = Raw_F_Ice(1); % Specify the inital value of the outlet flowrate
g_F_UG(1)    = Raw_F_UG(1);
g_L2(1)      = Raw_L(1);     % Specify the inital value of the level
for i = 1:51086            % From here onwards, make predictions
    e_F_outPT(i)   = randn*sqrt(w_F_outPT);
    g_F_outPT(i+1) = c_F_outPT + a_F_outPT*k_F_outPT(i) + e_F_outPT(i);
    e_F_Ice(i)     = randn*sqrt(w_F_Ice/100);
    g_F_Ice(i+1)   = c_F_Ice + a_F_Ice*k_F_Ice(i) + e_F_Ice(i);
    e_F_UG(i)      = randn*sqrt(w_F_UG);
    g_F_UG(i+1)    = c_F_UG + a_F_UG*k_F_UG(i) + e_F_UG(i);
    e_L(i)         = randn*sqrt(w_L);
    g_L2(i+1)      = a_L*k_L2(i) + C_L*(k_F_outPT(i) - k_F_Ice(i) - k_F_UG(i)) + e_L(i) + c_L;
end
x_K2 = [k_L2 ; k_F_outPT ; k_F_Ice; k_F_UG]; % State vector estimate/model prediction


g_F_Ice = g_F_Ice';
Size_g_F_Ice = size(g_F_Ice,1);
ipt_g = findchangepts(g_F_Ice,'Statistic','std','MinThreshold',100);
intvls_g = [1 ipt_g' length(g_F_Ice)];
for k1_g = 1:length(intvls_g)-2
    interval_means_g(k1_g) = mean(g_F_Ice(intvls_g(k1_g):intvls_g(k1_g+1)-1));
end
interval_means_g(k1_g+1) = mean(g_F_Ice(intvls_g(k1_g):intvls_g(end)));
xval_g = [];
yval_g = [];
for k1_g = 1:length(intvls_g)-1
    xval_g = intvls_g(k1_g):intvls_g(k1_g+1);
ynext_g = interval_means_g(k1_g)*ones(1,(length(xval_g)-1));
yval_g = [yval_g ynext_g];
end
yval_g(end+1) = 0;
g_F_Ice = yval_g;

figure(1)
plot(t,Raw_F_Ice);
hold on
plot(t,yval','LineWidth',2);
hold off
% INLET FLOWRATE
% Plot Actual Data vs AR Model

subplot(3,1,1)
plot(t, Raw_F_in, 'r', t, k_F_in', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('F_i_n_C_D (L/s)')
title('Inlet Flowrate AR Model');

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
% % Plot Actual Data vs AR Model
figure(2)
subplot(3,1,1)
plot(t, Raw_F_out, 'r', t, k_F_out', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('F_o_u_t_C_D (L/s)')
title('Outlet Flowrate AR Model');

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

% PRE-COOLING TOWER OUTLET FLOWRATE
% Plot Actual Data vs AR Model
figure(3)
subplot(3,1,1)
plot(t, Raw_F_outPT, 'r', t, k_F_outPT', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('F_o_u_t_P_T (L/s)')
title('Pre-Cooling Tower Outlet Flowrate AR Model');

% Plot Observations and Fitted Values
subplot(3,1,2)
resid_F_outPT = infer(EstMdl_F_outPT, Raw_F_outPT(estsample_F_outPT), ...
             'Y0', Raw_F_outPT(presample_F_outPT));        % Infer residuals from
                                                       % the estimated model
yhat_F_outPT = Raw_F_outPT(estsample_F_outPT) - resid_F_outPT; % Compute the fitted values
plot(t1, Raw_F_outPT(estsample_F_outPT),'r', t1, yhat_F_outPT,'b--','LineWidth',1)
legend('Observations', 'Fitted Values')

% Plot Residuals vs Fitted Values
subplot(3,1,3)
plot(yhat_F_outPT,resid_F_outPT,'b.')
ylabel('Residuals')
xlabel('Fitted Values')

% ICE PLANT FLOWRATE
% Plot Actual Data vs AR Model
figure(4)
subplot(3,1,1)
plot(t, Raw_F_Ice, 'r', t, k_F_Ice', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('F_I_c_e (L/s)')
title('Ice Plant Outlet Flowrate AR Model');

% Plot Observations and Fitted Values
subplot(3,1,2)
resid_F_Ice = infer(EstMdl_F_Ice, Raw_F_Ice(estsample_F_Ice), ...
             'Y0', Raw_F_Ice(presample_F_Ice));        % Infer residuals from
                                                       % the estimated model
yhat_F_Ice = Raw_F_Ice(estsample_F_Ice) - resid_F_Ice; % Compute the fitted values
plot(t1, Raw_F_Ice(estsample_F_Ice),'r', t1, yhat_F_Ice,'b--','LineWidth',1)
legend('Observations', 'Fitted Values')

% Plot Residuals vs Fitted Values
subplot(3,1,3)
plot(yhat_F_Ice,resid_F_Ice,'b.')
ylabel('Residuals')
xlabel('Fitted Values')

% UNDERGROUND FLOWRATE
% Plot Actual Data vs AR Model
figure(5)
subplot(3,1,1)
plot(t, Raw_F_UG, 'r', t, k_F_UG', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('F_U_G (L/s)')
title('UG Outlet Flowrate AR Model');

% Plot Observations and Fitted Values
subplot(3,1,2)
resid_F_UG = infer(EstMdl_F_UG, Raw_F_UG(estsample_F_UG), ...
             'Y0', Raw_F_UG(presample_F_UG));        % Infer residuals from
                                                       % the estimated model
yhat_F_UG = Raw_F_UG(estsample_F_UG) - resid_F_UG; % Compute the fitted values
plot(t1, Raw_F_UG(estsample_F_UG),'r', t1, yhat_F_UG,'b--','LineWidth',1)
legend('Observations', 'Fitted Values')

% Plot Residuals vs Fitted Values
subplot(3,1,3)
plot(yhat_F_UG,resid_F_UG,'b.')
ylabel('Residuals')
xlabel('Fitted Values')

% LEVEL FROM BALANCE OVER CHILL DAM
% Plot Actual Data vs AR Model
figure(6)
subplot(3,1,1)
plot(t, Raw_L, 'r', t, k_L1', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('L_C_D (m)')
title('Level AR Model from Balance over Chill Dam only');

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

% LEVEL FROM BALANCE OVER WHOLE SYSTEM
% Plot Actual Data vs AR Model
figure(7)
subplot(3,1,1)
plot(t, Raw_L, 'r', t, k_L2', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('L_C_D (m)')
title('Level AR Model from Balance over Whole System');

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
figure(8)
title('Level AR Model in %');
k_L_percent    = k_L1*100/p.height_CD;
L_Data_percent = Raw_L*100/p.height_CD;
plot(t, L_Data_percent, 'r', t, k_L_percent', 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (s)')
ylabel('L_C_D (%)')

save ArimaModelsCD.mat c_F_in c_F_out c_L w_F_in w_F_out w_L ...
    C_L a_F_in a_F_out a_L Raw_F_in ...
    Raw_F_out Raw_L c_F_outPT c_F_Ice c_F_UG a_F_outPT a_F_Ice ...
    a_F_UG w_F_outPT w_F_Ice w_F_UG Raw_F_outPT Raw_F_Ice Raw_F_UG...
    k_L2 k_F_outPT k_F_Ice k_F_UG k_L1 k_F_in k_F_out g_L1 g_F_out g_F_in...
    g_F_outPT g_F_Ice g_F_UG g_L2
