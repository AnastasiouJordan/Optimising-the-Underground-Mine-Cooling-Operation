function [p,u] = ARModels(t, u, p)

% ARIMA MODEL FOR T_RPj
Raw.T = u.T_RP(t);
Size_T = size(Raw.T,1);
t1 = t(1:51086,1);

% Create Model Template
Mdl_T = arima(1,0,0); 

presample_T = 1:Mdl_T.P;            % Contains 2 observations
estsample_T = (Mdl_T.P + 1):Size_T; % Contains the remaining observations

% Estimate the Model
% This is where we fit the ARIMA(1,1,0) Model to the estimation sample
for i = 1:5
    EstMdl_T(i) = estimate(Mdl_T, Raw.T(estsample_T,i), 'Y0', Raw.T(presample_T,i));
    p.a_T(i) = cell2mat(EstMdl_T(i).AR(1)); % AR coefficient, alpha
    p.w_T(i) = EstMdl_T(i).Variance;        % Variance which is rooted to obtain error
    p.c_T(i) = EstMdl_T(i).Constant;        % Constant or intercept
    k.T(1,i)  = Raw.T(1,i);  % Specify the inital value of the inlet flowrate
end

for j = 1:51086 % From here onwards, make predictions
    e.T(j,:) = randn.*sqrt(p.w_T);
    k.T(j+1,:) = p.c_T + p.a_T.*k.T(j,:) + e.T(j,:);
end

u.T_generated = griddedInterpolant(t, k.T, "previous");


% ARIMA MODEL FOR Q_refr
% Load the data
Raw.Q_refr  = (p.UA_RP .* u.T_refr(t)) ./ (p.C_p .* p.m_RPj);
Size_Q_refr= size(Raw.Q_refr,1);
t1 = t(1:51086,1);

% Create Model Template
dep = 0.998; % Dependence of the next L measurement on the previous one
Mdl_Q_refr = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant

%Mdl_Q_refr = arima(1,0,0); 

presample_Q_refr = 1:Mdl_Q_refr.P;                 % Contains 2 observations
estsample_Q_refr = (Mdl_Q_refr.P + 1):Size_Q_refr; % Contains the remaining observations

for i = 1:5
    EstMdl_Q_refr(i) = estimate(Mdl_Q_refr, Raw.Q_refr(estsample_Q_refr,i), 'Y0', Raw.Q_refr(presample_Q_refr,i));
    p.a_Q_refr(i) = cell2mat(EstMdl_Q_refr(i).AR(1)); % AR coefficient, alpha
    p.w_Q_refr(i) = EstMdl_Q_refr(i).Variance;        % Variance which is rooted to obtain error
    p.c_Q_refr(i) = EstMdl_Q_refr(i).Constant;        % Constant or intercept
    k.Q_refr(1,i)  = Raw.Q_refr(1,i);  % Specify the inital value of the inlet flowrate
end

for j = 1:51086 % From here onwards, make predictions
    e.Q_refr(j,:) = randn.*sqrt(p.w_Q_refr);
    k.Q_refr(j+1,:) = p.c_Q_refr + p.a_Q_refr.*k.Q_refr(j,:) + e.Q_refr(j,:);
end

u.Q_refr_generated = griddedInterpolant(t, k.Q_refr, "previous");


% for i = 1:5
%     figure(i)
%     plot(t, Raw.Q_refr(:,i), 'r', t, u.Q_refr_generated.Values(:,i)', 'b')
%     legend('Raw Data', 'AR Model')
%     xlabel('Time (s)')
%     ylabel('Q_refr')
% end


% % Model for T_amb
% y = u.T_amb(t);
% time_stamps = t;
% [xData, yData] = prepareCurveData( time_stamps, y );
% 
% % Set up fittype and options.
% ft = fittype( 'sin8' );
% opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
% opts.Display = 'Off';
% opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf];
% opts.StartPoint = [23.7648385274952 1.02493594252496e-06 0.0372165803905402 8.18264096413827 2.04987188504991e-06 1.66687793335002 4.31882657065534 3.48478220458485e-05 1.58008568665763 2.72281891504735 3.68976939308984e-05 1.27794503028682 2.25502717335074 3.27979501607986e-05 -1.10402109352573 2.72189470669459 4.09974377009982e-06 2.39142473064352 1.96830544550466 6.14961565514973e-06 1.25991284769313 1.51073542963859 7.17455159767469e-05 1.96470715835183];
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% % % Plot fit with data.
% % figure( 'Name', 'untitled fit 1' );
% % h = plot( fitresult, xData, yData );
% % legend( h, 'y vs. time_stamps', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % % Label axes
% % xlabel( 'time_stamps', 'Interpreter', 'none' );
% % ylabel( 'y', 'Interpreter', 'none' );
% % grid on
% 
% u.T_amb_generated = griddedInterpolant(t, fitresult(y), "previous");
% 


% ARIMA MODEL FOR Q_amb
% Load the data
Raw.Q_amb  = (p.UA_amb .* u.T_amb(t)) ./ (p.C_p .* p.m_RPj);
Size_Q_amb= size(Raw.Q_amb,1);
t1 = t(1:51086,1);

% Create Model Template
%Mdl_Q_amb = arima(1,0,0); 
dep = 0.9997; % Dependence of the next L measurement on the previous one
Mdl_Q_amb = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant

presample_Q_amb = 1:Mdl_Q_amb.P;                 % Contains 2 observations
estsample_Q_amb = (Mdl_Q_amb.P + 1):Size_Q_amb; % Contains the remaining observations

for i = 1:5
    EstMdl_Q_amb(i) = estimate(Mdl_Q_amb, Raw.Q_amb(estsample_Q_amb,i), 'Y0', Raw.Q_amb(presample_Q_amb,i));
    p.a_Q_amb(i) = cell2mat(EstMdl_Q_amb(i).AR(1)); % AR coefficient, alpha
    p.w_Q_amb(i) = EstMdl_Q_amb(i).Variance;        % Variance which is rooted to obtain error
    p.c_Q_amb(i) = EstMdl_Q_amb(i).Constant;        % Constant or intercept
    k.Q_amb(1,i)  = Raw.Q_amb(1,i);  % Specify the inital value of the inlet flowrate
end

for j = 1:51086 % From here onwards, make predictions
    e.Q_amb(j,:) = randn.*sqrt(p.w_Q_amb);
    k.Q_amb(j+1,:) = p.c_Q_amb + p.a_Q_amb.*k.Q_amb(j,:) + e.Q_amb(j,:);
end

u.Q_amb_generated = griddedInterpolant(t, k.Q_amb, "previous");


% for i = 1:5
%     figure(i)
%     plot(t, Raw.Q_amb(:,i), 'r', t, u.Q_amb_generated.Values(:,i)', 'b')
%     legend('Raw Data', 'AR Model')
%     xlabel('Time (s)')
%     ylabel('Q_amb')
% end

% ARIMA MODEL FOR T_inRPj
Raw.T_in = u.T_inRP(t);
Size_T_in = size(Raw.T_in,1);
t1 = t(1:51086,1);

% Create Model Template
Mdl_T_in = arima(1,0,0); 

presample_T_in = 1:Mdl_T_in.P;            % Contains 2 observations
estsample_T_in = (Mdl_T_in.P + 1):Size_T_in; % Contains the remaining observations

for i = 1:5
    EstMdl_T_in(i) = estimate(Mdl_T_in, Raw.T_in(estsample_T_in,i), 'Y0', Raw.T_in(presample_T_in,i));
    p.a_T_in(i) = cell2mat(EstMdl_T_in(i).AR(1)); % AR coefficient, alpha
    p.w_T_in(i) = EstMdl_T_in(i).Variance;        % Variance which is rooted to obtain error
    p.c_T_in(i) = EstMdl_T_in(i).Constant;        % Constant or intercept
    % AR Model Predictions
    k.T_in(1,i)  = Raw.T_in(1,i);  % Specify the inital value of the inlet flowrate
end

for j = 1:51086 % From here onwards, make predictions
    e.T_in(j,:) = randn.*sqrt(p.w_T_in);
    k.T_in(j+1,:) = p.c_T_in + p.a_T_in.*k.T_in(j,:) + e.T_in(j,:);
end

u.T_in_generated = griddedInterpolant(t, k.T_in, "previous");

% plot(t, Raw.T_in(:,i), 'r', t, u.T_in_generated.Values(:,i)', 'b')
%     legend('Raw Data', 'AR Model')
%     xlabel('Time (s)')
%     ylabel('T_in')

end