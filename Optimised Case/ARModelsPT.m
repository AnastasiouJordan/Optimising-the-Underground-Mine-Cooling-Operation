function [p,u] = ARModelsPT(t, u, p)

% ARIMA Models

% ARIMA MODEL FOR F_inPT
% Load the data
Raw.F_in  = u.F_outSD(t) + u.F_AntiSurge(t);
Size_F_in = size(Raw.F_in,1);
t1 = t(1:51086,1);

% Create Model Template
%Mdl_F_in = arima(1,0,0); % First-order Auto-Regressive Model plus constant
dep = 0.97; % Dependence of the next L measurement on the previous one
Mdl_F_in = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_F_in = 1:Mdl_F_in.P;               % Contains 2 observations
estsample_F_in = (Mdl_F_in.P + 1):Size_F_in; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_F_in = estimate(Mdl_F_in, Raw.F_in(estsample_F_in), 'Y0', Raw.F_in(presample_F_in));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_F_in_PT = cell2mat(EstMdl_F_in.AR(1)); % AR coefficient, alpha
p.w_F_in_PT = EstMdl_F_in.Variance;        % Variance which is rooted to obtain error
p.c_F_in_PT = EstMdl_F_in.Constant;        % Constant or intercept




% ARIMA MODEL FOR F_Antisurge
% Load the data
Raw.F_Antisurge  = u.F_AntiSurge(t);
Size_F_Antisurge = size(Raw.F_Antisurge,1);
t1 = t(1:51086,1);
% Create Model Template
%Mdl_F_Antisurge = arima(1,0,0); % First-order Auto-Regressive Model plus constant
dep = 0.999; % Dependence of the next L measurement on the previous one
Mdl_F_Antisurge = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant

% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_F_Antisurge = 1:Mdl_F_Antisurge.P;               % Contains 2 observations
estsample_F_Antisurge = (Mdl_F_Antisurge.P + 1):Size_F_Antisurge; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_F_Antisurge = estimate(Mdl_F_Antisurge, Raw.F_Antisurge(estsample_F_Antisurge), 'Y0', Raw.F_Antisurge(presample_F_Antisurge));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_F_Antisurge = cell2mat(EstMdl_F_Antisurge.AR(1)); % AR coefficient, alpha
p.w_F_Antisurge = EstMdl_F_Antisurge.Variance;        % Variance which is rooted to obtain error
p.c_F_Antisurge = EstMdl_F_Antisurge.Constant;        % Constant or intercept



% ARIMA MODEL FOR T_Antisurge
% Load the data
Raw.T_Antisurge  = u.T_AntiSurge(t);
Size_T_Antisurge = size(Raw.T_Antisurge,1);
t1 = t(1:51086,1);

% Create Model Template
%Mdl_T_Antisurge = arima(1,0,0); % First-order Auto-Regressive Model plus constant
dep = 0.9995; % Dependence of the next L measurement on the previous one
Mdl_T_Antisurge = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant

% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_T_Antisurge = 1:Mdl_T_Antisurge.P;               % Contains 2 observations
estsample_T_Antisurge = (Mdl_T_Antisurge.P + 1):Size_T_Antisurge; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_T_Antisurge = estimate(Mdl_T_Antisurge, Raw.T_Antisurge(estsample_T_Antisurge), 'Y0', Raw.T_Antisurge(presample_T_Antisurge));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_T_Antisurge = cell2mat(EstMdl_T_Antisurge.AR(1)); % AR coefficient, alpha
p.w_T_Antisurge = EstMdl_T_Antisurge.Variance;        % Variance which is rooted to obtain error
p.c_T_Antisurge = EstMdl_T_Antisurge.Constant;        % Constant or intercept



% ARIMA MODEL FOR F_outPT

% Load the data
Raw.F_out = u.F_outPT(t);
Size_F_out = size(Raw.F_out,1);

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
EstMdl_F_out = estimate(Mdl_F_out, Raw.F_out(estsample_F_out), 'Y0', Raw.F_out(presample_F_out));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_F_out_PT = cell2mat(EstMdl_F_out.AR(1)); % AR coefficient, alpha
p.w_F_out_PT = EstMdl_F_out.Variance;        % Variance which is rooted to obtain error
p.c_F_out_PT = EstMdl_F_out.Constant;        % Constant or intercept


% ARIMA MODEL FOR L_PT

% Load the data
Raw_L = (u.L_PTA(t) + u.L_PTB(t))./2; % Raw level measurements as percentages
Raw.L = Raw_L/100*p.height_PT;        % Raw level measurements as heighs in meters
Size_L = size(Raw.L,1);

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
EstMdl_L = estimate(Mdl_L, Raw.L(estsample_L), 'Y0', Raw.L(presample_L));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_L_PT    = cell2mat(EstMdl_L.AR(1)); % AR coefficient, alpha
p.w_L_PT    = EstMdl_L.Variance;        % Variance which is rooted to obtain error
p.c_L_PT    = EstMdl_L.Constant;        % Constant or intercept
deltaT = t(2) - t(1);                % s, Delta time in seconds
p.C_L_PT    = (deltaT/p.area_PT)*0.001;   % m, Coefficient to the flowrates.
                                     % This calculation takes delta time and
                                     % divides it by the area of the Chill Dam
                                     % in m2. This is then multiplied by 0.001
                                     % to adjust for the flowrate which is in
                                     % L/s. This gives a constant in meters,
                                     % equivalent to the units in which level
                                     % is measured.





% ARIMA MODEL FOR T_PT

% Load the data
Raw.T_PT = u.T_PT(t); % Raw level measurements as percentages
Size_T_PT = size(Raw.T_PT,1);

% Create Model Template
dep = 0.997; % Dependence of the next L measurement on the previous one
%Mdl_L_A = arima(1,0,0); % This gives an alpha value of 0.9999 and poor fit
Mdl_T_PT = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_T_PT = 1:Mdl_T_PT.P;            % Contains 2 observations
estsample_T_PT = (Mdl_T_PT.P + 1):Size_T_PT; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_T_PT = estimate(Mdl_T_PT, Raw.T_PT(estsample_T_PT), 'Y0', Raw.T_PT(presample_T_PT));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_T_PT    = cell2mat(EstMdl_T_PT.AR(1)); % AR coefficient, alpha
p.w_T_PT    = EstMdl_T_PT.Variance;        % Variance which is rooted to obtain error
p.c_T_PT    = EstMdl_T_PT.Constant;        % Constant or intercept



% ARIMA MODEL FOR H_Env

% Load the data
Raw.H_Env = u.H_Env(t); % Raw level measurements as percentages
Size_H_Env = size(Raw.H_Env,1);

% Create Model Template
dep = 0.9998; % Dependence of the next L measurement on the previous one
%Mdl_L_A = arima(1,0,0); % This gives an alpha value of 0.9999 and poor fit
Mdl_H_Env = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_H_Env = 1:Mdl_H_Env.P;            % Contains 2 observations
estsample_H_Env = (Mdl_H_Env.P + 1):Size_H_Env; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_H_Env = estimate(Mdl_H_Env, Raw.H_Env(estsample_H_Env), 'Y0', Raw.H_Env(presample_H_Env));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_H_Env    = cell2mat(EstMdl_H_Env.AR(1)); % AR coefficient, alpha
p.w_H_Env    = EstMdl_H_Env.Variance;        % Variance which is rooted to obtain error
p.c_H_Env    = EstMdl_H_Env.Constant;        % Constant or intercept



% AR Model Predictions
g.F_in(1)  = Raw.F_in(1);  % Specify the inital value of the inlet flowrate
g.F_out(1) = Raw.F_out(1); % Specify the inital value of the outlet flowrate
g.L(1)     = Raw.L(1);     % Specify the inital value of the level
g.F_Antisurge(1) = Raw.F_Antisurge(1);
g.T_Antisurge(1) = Raw.T_Antisurge(1);
g.T_PT(1) = Raw.T_PT(1);
g.H_Env(1) = Raw.H_Env(1);
for i = 1:51086           
    e.F_in(i)    = randn*sqrt(p.w_F_in_PT);
    g.F_in(i+1)  = p.c_F_in_PT + p.a_F_in_PT*g.F_in(i) + e.F_in(i);    % Generated data used to test the model in simulation
    e.F_Antisurge(i)   = randn*sqrt(p.w_F_Antisurge);
    g.F_Antisurge(i+1) = p.c_F_Antisurge + p.a_F_Antisurge*g.F_Antisurge(i) + e.F_Antisurge(i);    % Generated data used to test the model in simulation
    e.T_Antisurge(i)   = randn*sqrt(p.w_T_Antisurge);
    g.T_Antisurge(i+1) = p.c_T_Antisurge + p.a_T_Antisurge*g.T_Antisurge(i) + e.T_Antisurge(i);    % Generated data used to test the model in simulation
    e.F_out(i)   = randn*sqrt(p.w_F_out_PT);
    g.F_out(i+1) = p.c_F_out_PT + p.a_F_out_PT*g.F_out(i) + e.F_out(i); % Generated data used to test the model in simulation
    e.L(i)       = randn*sqrt(p.w_L_PT);
    g.L(i+1)     = p.a_L_PT*g.L(i) + p.C_L_PT*(g.F_in(i) - g.F_out(i)) + e.L(i) + p.c_L_PT; % Generated data used to test the model in simulation
    e.T_PT(i)   = randn*sqrt(p.w_T_PT);
    g.T_PT(i+1) = p.c_T_PT + p.a_T_PT*g.T_PT(i) + e.T_PT(i);    % Generated data used to test the model in simulation
    e.H_Env(i)   = randn*sqrt(p.w_H_Env);
    g.H_Env(i+1) = p.c_H_Env + p.a_H_Env*g.H_Env(i) + e.H_Env(i);    % Generated data used to test the model in simulation
end
% Where in the above loop:
% e(i) is the error obtained from the variance calculated by the model for each variable.
% c is the constant or intercept calculated by the model for each variable.
% a is the 'alpha' co-efficient for the dependence of the next value on the previous one for each variable.
% k(i) is the prediction for each variable.
% C is the constant used only in the level predictions, derived from the relationship between level and flowrates.


% g_F_in = g.F_in';
% Size_g_F_in = size(g_F_in,1);
% ipt_g = findchangepts(g_F_in,'Statistic','std','MinThreshold',150);
% intvls_g = [1 ipt_g' length(g_F_in)];
% for k1_g = 1:length(intvls_g)-2
%     interval_means_g(k1_g) = mean(g_F_in(intvls_g(k1_g):intvls_g(k1_g+1)-1));
% end
% interval_means_g(k1_g+1) = mean(g_F_in(intvls_g(k1_g):intvls_g(end)));
% xval_g = [];
% yval_g = [];
% for k1_g = 1:length(intvls_g)-1
%     xval_g = intvls_g(k1_g):intvls_g(k1_g+1);
% ynext_g = interval_means_g(k1_g)*ones(1,(length(xval_g)-1));
% yval_g = [yval_g ynext_g];
% end
% yval_g(end+1) = 0;
% g_F_in = yval_g;


g_F_Antisurge = g.F_Antisurge';
yval_g = [];
Size_g_Antisurge = size(g_F_Antisurge,1);
ipt_g = findchangepts(g_F_Antisurge,'Statistic','std','MinThreshold',600);
intvls_g = [1 ipt_g' length(g_F_Antisurge)];
for k1_g = 1:length(intvls_g)-2
    interval_means_g(k1_g) = mean(g_F_Antisurge(intvls_g(k1_g):intvls_g(k1_g+1)-1));
end
interval_means_g(k1_g+1) = mean(g_F_Antisurge(intvls_g(k1_g):intvls_g(end)));
xval_g = [];
yval_g = [];
for k1_g = 1:length(intvls_g)-1
    xval_g = intvls_g(k1_g):intvls_g(k1_g+1);
ynext_g = interval_means_g(k1_g)*ones(1,(length(xval_g)-1));
yval_g = [yval_g ynext_g];
end
yval_g(end+1) = 0;
g.F_Antisurge = 1000*yval_g;
g.F_Antisurge = max(0, g.F_Antisurge);
g.F_Antisurge = min(50, g.F_Antisurge);

u.F_in_generatedPT = griddedInterpolant(t, g.F_in', "previous");
u.F_out_generatedPT = griddedInterpolant(t, g.F_out');
u.F_AntiSurge_generated = griddedInterpolant(t, g.F_Antisurge');
u.T_AntiSurge_generated = griddedInterpolant(t, g.T_Antisurge');
u.T_PT_generated = griddedInterpolant(t, g.T_PT');


g.H_Env = max(0, g.H_Env);
g.H_Env = min(100, g.H_Env);
u.H_Env_generated = griddedInterpolant(t, g.H_Env');


days = 86400;

% % Model for T_amb
y = u.T_amb(t);
time_stamps = t;
[xData, yData] = prepareCurveData( time_stamps, y );

% Set up fittype and options.
ft = fittype( 'sin8' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf];
opts.StartPoint = [23.7648385274952 1.02493594252496e-06 0.0372165803905402 8.18264096413827 2.04987188504991e-06 1.66687793335002 4.31882657065534 3.48478220458485e-05 1.58008568665763 2.72281891504735 3.68976939308984e-05 1.27794503028682 2.25502717335074 3.27979501607986e-05 -1.10402109352573 2.72189470669459 4.09974377009982e-06 2.39142473064352 1.96830544550466 6.14961565514973e-06 1.25991284769313 1.51073542963859 7.17455159767469e-05 1.96470715835183];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
u.T_amb_generated = griddedInterpolant(t, fitresult(t), "previous");

Raw.T_out = u.T_outSD(t);
Size_T_out = size(Raw.T_out,1);
% Create Model Template
Mdl_T_out = arima(1,0,0); % First-order Auto-Regressive Model plus constant
presample_T_out = 1:Mdl_T_out.P;                % Contains 2 observations
estsample_T_out = (Mdl_T_out.P + 1):Size_T_out; % Contains the remaining observations
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_T_out = estimate(Mdl_T_out, Raw.T_out(estsample_T_out), 'Y0', Raw.T_out(presample_T_out));


% figure(1)
% plot(t/days, Raw.T_Antisurge, 'r', t/days, u.T_AntiSurge_generated(t)', 'b')
% legend('Raw Data', 'AR Model')
% xlabel('Time (days)')
% ylabel('T_A_n_t_i_S_u_r_g_e (^oC)')
% ax = gca;
% ax.FontSize = p.font_size;
% 
% figure(2)
% plot(t/days, Raw.F_Antisurge, 'r', t/days, u.F_AntiSurge_generated(t)', 'b')
% legend('Raw Data', 'AR Model')
% xlabel('Time (days)')
% ylabel('F_A_n_t_i_S_u_r_g_e (L/s)')
% ax = gca;
% ax.FontSize = p.font_size;
% 
% figure(3)
% plot(t/days, Raw.T_PT, 'r', t/days, u.T_PT_generated(t)', 'b')
% legend('Raw Data', 'AR Model')
% xlabel('Time (days)')
% ylabel('T_PT')
% ax = gca;
% ax.FontSize = p.font_size;
% 
% figure(4)
% plot(t/days, Raw.H_Env, 'r', t/days, u.H_Env_generated(t)', 'b')
% legend('Raw Data', 'AR Model')
% xlabel('Time (days)')
% ylabel('H_E_n_v (%)')
% ax = gca;
% ax.FontSize = p.font_size;
% 
% figure(5)
% plot(t/days, Raw.T_out, 'r', t/days, u.T_out_generatedSD(t)', 'b')
% legend('Raw Data', 'AR Model')
% xlabel('Time (days)')
% ylabel('T_o_u_t_S_D (^oC)')
% ax = gca;
% ax.FontSize = p.font_size;

%% Plotting Residuals

figure(1)
%T_Antisurge
subplot(2,3,1)
title('Anti-surge Temperature');
hold on
plot(t/86400, Raw.T_Antisurge, 'r', t/86400, u.T_AntiSurge_generated(t), 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (days)')
ylabel('T_A_n_t_i_S_u_r_g_e (^oC)')
ax = gca;
ax.FontSize = p.font_size;
hold off

% Plot Observations and Fitted Values
resid_T_Antisurge = infer(EstMdl_T_Antisurge, Raw.T_Antisurge(estsample_T_Antisurge), ...
             'Y0', Raw.T_Antisurge(presample_T_Antisurge));      % Infer residuals from
                                                   % the estimated model
yhat_T_Antisurge = Raw.T_Antisurge(estsample_T_Antisurge) - resid_T_Antisurge; % Compute the fitted values


% Plot Residuals vs Fitted Values
subplot(2,3,4)
plot(yhat_T_Antisurge,resid_T_Antisurge,'b.')
ylabel('Residuals')
xlabel('Fitted Values')
ax = gca;
ax.FontSize = p.font_size;

%T_outSD
subplot(2,3,2)
title('Inlet Temperature');
hold on
plot(t/86400, Raw.T_out, 'r', t/86400, u.T_out_generatedSD(t), 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (days)')
ylabel('T_o_u_t_S_D (^oC)')
ax = gca;
ax.FontSize = p.font_size;
hold off

% Plot Observations and Fitted Values
resid_T_out = infer(EstMdl_T_out, Raw.T_out(estsample_T_out), ...
             'Y0', Raw.T_out(presample_T_out));      % Infer residuals from
                                                   % the estimated model
yhat_T_out = Raw.T_out(estsample_T_out) - resid_T_out; % Compute the fitted values


% Plot Residuals vs Fitted Values
subplot(2,3,5)
plot(yhat_T_out,resid_T_out,'b.')
ylabel('Residuals')
xlabel('Fitted Values')
ax = gca;
ax.FontSize = p.font_size;

%H_Env
subplot(2,3,3)
title('Environmental Humidity');
hold on
plot(t/86400, Raw.H_Env, 'r', t/86400, u.H_Env_generated(t), 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (days)')
ylabel('H_E_n_v (%)')
ax = gca;
ax.FontSize = p.font_size;
hold off

% Plot Observations and Fitted Values
resid_H_Env = infer(EstMdl_H_Env, Raw.H_Env(estsample_H_Env), ...
             'Y0', Raw.H_Env(presample_H_Env));      % Infer residuals from
                                                   % the estimated model
yhat_H_Env = Raw.H_Env(estsample_H_Env) - resid_H_Env; % Compute the fitted values


% Plot Residuals vs Fitted Values
subplot(2,3,6)
plot(yhat_H_Env,resid_H_Env,'b.')
ylabel('Residuals')
xlabel('Fitted Values')
ax = gca;
ax.FontSize = p.font_size;

%% T amb and F Antisurge
figure(2)
subplot(2,1,1)
title('Anti-surge Flowrate');
hold on
plot(t/86400, Raw.F_Antisurge, 'r', t/86400, u.F_AntiSurge_generated(t), 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (days)')
ylabel('F_A_n_t_i_S_u_r_g_e (L/s)')
ax = gca;
ax.FontSize = p.font_size;
hold off

Raw.T_amb = u.T_amb(t);
subplot(2,1,2)
title('Ambient Temperature');
hold on
plot(t/86400, Raw.T_amb, 'r', t/86400, u.T_amb_generated(t), 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (days)')
ylabel('T_a_m_b (^oC)')
ax = gca;
ax.FontSize = p.font_size;
hold off
end
