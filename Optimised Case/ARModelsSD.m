function [p,u] = ARModelsSD(t, u, p)

% ARIMA Models

% ARIMA MODEL FOR F_inSD
% Load the data
Raw.F_in  = u.F_inSD(t);
Size_F_in = size(Raw.F_in,1);
t1 = t(1:51086,1);

% Create Model Template
dep = 0.98; % Dependence of the next L measurement on the previous one
Mdl_F_in = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant

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
EstMdl_F_in = estimate(Mdl_F_in, Raw.F_in(estsample_F_in), 'Y0', Raw.F_in(presample_F_in));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_F_in_SD = cell2mat(EstMdl_F_in.AR(1)); % AR coefficient, alpha
p.w_F_in_SD = EstMdl_F_in.Variance;        % Variance which is rooted to obtain error
p.c_F_in_SD = EstMdl_F_in.Constant;        % Constant or intercept


% ARIMA MODEL FOR F_outSD

% Load the data
Raw.F_out = u.F_outSD(t);
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
p.a_F_out_SD = cell2mat(EstMdl_F_out.AR(1)); % AR coefficient, alpha
p.w_F_out_SD = EstMdl_F_out.Variance;        % Variance which is rooted to obtain error
p.c_F_out_SD = EstMdl_F_out.Constant;        % Constant or intercept


% ARIMA MODEL FOR L_SD

% Load the data
Raw.L = u.L_SD(t); % Raw level measurements as percentages
Raw.L = Raw.L/100*p.height_SD; % Raw level measurements as heights in meters
Size_L = size(Raw.L,1);

% Create Model Template
dep = 0.99; % Dependence of the next L measurement on the previous one
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
p.a_L_SD    = cell2mat(EstMdl_L.AR(1)); % AR coefficient, alpha
p.w_L_SD    = EstMdl_L.Variance;        % Variance which is rooted to obtain error
p.c_L_SD    = EstMdl_L.Constant;        % Constant or intercept
deltaT = t(2) - t(1);                % s, Delta time in seconds
p.C_L_SD    = (deltaT/p.area_SD)*0.001;   % m, Coefficient to the flowrates.
                                     % This calculation takes delta time and
                                     % divides it by the area of the Chill Dam
                                     % in m2. This is then multiplied by 0.001
                                     % to adjust for the flowrate which is in
                                     % L/s. This gives a constant in meters,
                                     % equivalent to the units in which level
                                     % is measured.



% ARIMA MODEL FOR T_outSD

% Load the data
Raw.T_out = u.T_outSD(t);
Size_T_out = size(Raw.T_out,1);

% Create Model Template
Mdl_T_out = arima(1,0,0); % First-order Auto-Regressive Model plus constant


% Partition Sample
% Creating vectors of indices that partition the sample into a 'presample'
% and an 'estimation sample' period
presample_T_out = 1:Mdl_T_out.P;                % Contains 2 observations
estsample_T_out = (Mdl_T_out.P + 1):Size_T_out; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
EstMdl_T_out = estimate(Mdl_T_out, Raw.T_out(estsample_T_out), 'Y0', Raw.T_out(presample_T_out));

% Save Equation Constants
% The AR model is based on the following equation:
% k,i+1 = c + a*k,i + w
p.a_T_out_SD = cell2mat(EstMdl_T_out.AR(1)); % AR coefficient, alpha
p.w_T_out_SD = EstMdl_T_out.Variance;        % Variance which is rooted to obtain error
p.c_T_out_SD = EstMdl_T_out.Constant;        % Constant or intercept




% AR Model Predictions
k.F_in(1)  = Raw.F_in(1);  % Specify the inital value of the inlet flowrate
k.F_out(1) = Raw.F_out(1); % Specify the inital value of the outlet flowrate
k.L(1)     = Raw.L(1);     % Specify the inital value of the level
k.T_out(1) = Raw.T_out(1);
for i = 1:51086            % From here onwards, make predictions
    e.F_in(i)    = randn*sqrt(p.w_F_in_SD);
    k.F_in(i+1)  = p.c_F_in_SD + p.a_F_in_SD*k.F_in(i) + e.F_in(i);
    e.F_out(i)   = randn*sqrt(p.w_F_out_SD);
    k.F_out(i+1) = p.c_F_out_SD + p.a_F_out_SD*k.F_out(i) + e.F_out(i);
    e.T_out(i)   = randn*sqrt(p.w_T_out_SD);
    k.T_out(i+1) = p.c_T_out_SD + p.a_T_out_SD*k.T_out(i) + e.T_out(i);
    e.L(i)       = randn*sqrt(p.w_L_SD);
    k.L(i+1)     = p.a_L_SD*k.L(i) + p.C_L_SD*(k.F_in(i) - k.F_out(i)) + e.L(i) + p.c_L_SD;
end

g.F_in(1)  = Raw.F_in(1);  % Specify the inital value of the inlet flowrate
g.F_out(1) = Raw.F_out(1); % Specify the inital value of the outlet flowrate
g.L(1)     = Raw.L(1);     % Specify the inital value of the level
g.T_out(1) = Raw.T_out(1);
for i = 1:51086           
    e.F_in(i)    = randn*sqrt(p.w_F_in_SD);
    g.F_in(i+1)  = p.c_F_in_SD + p.a_F_in_SD*g.F_in(i) + e.F_in(i);    % Generated data used to test the model in simulation
    e.F_out(i)   = randn*sqrt(p.w_F_out_SD);
    g.F_out(i+1) = p.c_F_out_SD + p.a_F_out_SD*g.F_out(i) + e.F_out(i); % Generated data used to test the model in simulation
    e.T_out(i)   = randn*sqrt(p.w_T_out_SD);
    g.T_out(i+1) = p.c_T_out_SD + p.a_T_out_SD*g.T_out(i) + e.T_out(i);
    e.L(i)       = randn*sqrt(p.w_L_SD);
    g.L(i+1)     = p.a_L_SD*g.L(i) + p.C_L_SD*(g.F_in(i) - g.F_out(i)) + e.L(i) + p.c_L_SD; % Generated data used to test the model in simulation
end
% Where in the above loop:
% e(i) is the error obtained from the variance calculated by the model for each variable.
% c is the constant or intercept calculated by the model for each variable.
% a is the 'alpha' co-efficient for the dependence of the next value on the previous one for each variable.
% k(i) is the prediction for each variable.
% C is the constant used only in the level predictions, derived from the relationship between level and flowrates.




u.F_in_generatedSD = griddedInterpolant(t, g.F_in');
u.F_out_generatedSD = griddedInterpolant(t, g.F_out');
u.T_out_generatedSD = griddedInterpolant(t, g.T_out');

% plot(t/86400, Raw.F_in, 'r', t/86400, u.F_in_generatedSD(t)', 'b')
%     legend('Raw Data', 'AR Model')
%     xlabel('Time (days)')
%     ylabel('F_i_n_S_D (L/s)')
% ax = gca;
% ax.FontSize = p.font_size;

%% Plotting Residuals

%F_inSD
subplot(2,1,1)
title('Inlet Flowrate');
hold on
plot(t/86400, Raw.F_in, 'r', t/86400, u.F_in_generatedSD(t), 'b')
legend('Raw Data', 'AR Model')
xlabel('Time (days)')
ylabel('F_i_n_S_D (L/s)')
ax = gca;
ax.FontSize = p.font_size;
hold off

% Plot Observations and Fitted Values
resid_F_in = infer(EstMdl_F_in, Raw.F_in(estsample_F_in), ...
             'Y0', Raw.F_in(presample_F_in));      % Infer residuals from
                                                   % the estimated model
yhat_F_in = Raw.F_in(estsample_F_in) - resid_F_in; % Compute the fitted values
% plot(t1, Raw_F_in(estsample_F_in),'r', t1, yhat_F_in,'b--','LineWidth',1)
% legend('Observations', 'Fitted Values')

% Plot Residuals vs Fitted Values
subplot(2,1,2)
plot(yhat_F_in,resid_F_in,'b.')
ylabel('Residuals')
xlabel('Fitted Values')
ax = gca;
ax.FontSize = p.font_size;

end
