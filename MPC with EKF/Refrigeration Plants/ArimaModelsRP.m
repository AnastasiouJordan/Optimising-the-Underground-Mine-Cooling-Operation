clc
clear

load SavedInterpolantsRP.mat

% ARIMA Models

% % Define Dimensions
% p.height_PT = 3;     % m, Height of the Pre-cooling Tower basins
% p.length_PT = 14.8;  % m, Length of the Pre-cooling Tower basins
% p.width_PT  = 10.45*2;  % m, Width of the Pre-cooling Tower basins
% p.area_PT   = p.length_PT*p.width_PT; % m2, Area of the Pre-cooling Tower basins
% p.volume_PT = p.area_PT*p.height_PT;  % m3, Volume of the Pre-cooling Tower basins


% ARIMA MODEL FOR F_inRP
% Load the data
Raw_F_in  = u.F_inRPtot(t);
Size_F_in = size(Raw_F_in,1);
t1 = t(1:51086,1);


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


% AR Model Predictions
k_F_in(1)  = Raw_F_in(1);  % Specify the inital value of the inlet flowrate
for i = 1:51086            % From here onwards, make predictions
    e_F_in(i)    = randn*sqrt(w_F_in);
    k_F_in(i+1)  = c_F_in + a_F_in*k_F_in(i) + e_F_in(i);
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
ylabel('F_i_n_R_P_t_o_t (L/s)')

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


u.F_inRPtot_generated = griddedInterpolant(t, k_F_in');

save ArimaModelsRP.mat c_F_in w_F_in a_F_in k_F_in Raw_F_in u
