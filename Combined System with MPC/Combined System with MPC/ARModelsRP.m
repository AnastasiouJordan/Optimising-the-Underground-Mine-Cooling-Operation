function [p,u] = ARModelsRP(t, u, p)

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
    p.a_T_RP(i) = cell2mat(EstMdl_T(i).AR(1)); % AR coefficient, alpha
    p.w_T_RP(i) = EstMdl_T(i).Variance;        % Variance which is rooted to obtain error
    p.c_T_RP(i) = EstMdl_T(i).Constant;        % Constant or intercept
    k.T(1,i)  = Raw.T(1,i);  % Specify the inital value of the inlet flowrate
end

for j = 1:51086 % From here onwards, make predictions
    e.T(j,:) = randn.*sqrt(p.w_T_RP);
    k.T(j+1,:) = p.c_T_RP + p.a_T_RP.*k.T(j,:) + e.T(j,:);
end

u.T_generatedRP = griddedInterpolant(t, k.T, "previous");


% ARIMA MODEL FOR Q_refr
% Load the data
Raw.Q_refr  = min(1,(p.UA_RP .* u.T_refr(t)) ./ (p.C_p .* p.m_RPj));
Size_Q_refr= size(Raw.Q_refr,1);
t1 = t(1:51086,1);

% Create Model Template
dep = 0.95; % Dependence of the next L measurement on the previous one
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

k.Q_refr(:,2) = k.Q_refr(:,5);


u.Q_refr_generated = griddedInterpolant(t, k.Q_refr, "previous");



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

% Plot fit with data.
% figure(6);
% plot(fitresult, xData, yData);
% legend('Raw Data', 'Sin Function Fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'Time (days)');
% ylabel( 'T_a_m_b (^oC)');
% ax = gca;
% ax.FontSize = p.font_size;

u.T_amb_generated = griddedInterpolant(t, fitresult(t), "previous");



% ARIMA MODEL FOR Q_amb
% Load the data
Raw.Q_amb  = min(-0.05,(p.UA_amb .* u.T_amb(t)) ./ (p.C_p .* p.m_RPj));
Size_Q_amb= size(Raw.Q_amb,1);
t1 = t(1:51086,1);

% Create Model Template
%Mdl_Q_amb = arima(1,0,0); 
dep = 0.9995; % Dependence of the next L measurement on the previous one
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
k.Q_refr(:,2) = k.Q_refr(:,3);
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
    p.a_T_in_RP(i) = cell2mat(EstMdl_T_in(i).AR(1)); % AR coefficient, alpha
    p.w_T_in_RP(i) = EstMdl_T_in(i).Variance;        % Variance which is rooted to obtain error
    p.c_T_in_RP(i) = EstMdl_T_in(i).Constant;        % Constant or intercept
    % AR Model Predictions
    k.T_in(1,i)  = Raw.T_in(1,i);  % Specify the inital value of the inlet flowrate
end

for j = 1:51086 % From here onwards, make predictions
    e.T_in(j,:) = randn.*sqrt(p.w_T_in_RP);
    k.T_in(j+1,:) = p.c_T_in_RP + p.a_T_in_RP.*k.T_in(j,:) + e.T_in(j,:);
end

u.T_in_generatedRP = griddedInterpolant(t, k.T_in, "previous");


% for i = 1:5
%     figure(i)
%     plot(t/86400, Raw.Q_refr(:,i), 'r', t/86400, u.Q_refr_generated.Values(:,i)', 'b')
%     legend('Raw Data', 'AR Model')
%     xlabel('Time (days)')
%     ylabel('Q_r_e_f_r (^oC/s)')
%     ax = gca;
%     ax.FontSize = p.font_size;
% end


% figure(6)
% Raw.T_amb = u.T_amb(t);
% plot(t/86400, Raw.T_amb, 'r', t/86400, u.T_amb_generated(t), 'b')
% legend('Raw Data', 'Sin Function Fit')
% xlabel('Time (days)')
% ylabel('T_a_m_b (^oC)')
% ax = gca;
% ax.FontSize = p.font_size;

% for i = 1:5
%     figure(i)
%     plot(t/86400, Raw.Q_amb(:,i), 'r', t/86400, u.Q_amb_generated.Values(:,i)', 'b')
%     legend('Raw Data', 'AR Model')
%     xlabel('Time (days)')
%     ylabel('Q_a_m_b (^oC/s)')
%     ax = gca;
%     ax.FontSize = p.font_size;
% end





% Raw.T_amb_test = u.T_amb(t); % Raw level measurements as percentages
% Size_T_amb_test = size(Raw.T_amb_test,1);
% Mdl_T_amb_test = arima('Constant',1,'D',1,'Seasonality',12,...
% 	'MALags',1,'SMALags',12); % This gives an alpha value of 0.9999 and poor fit
% presample_T_amb_test = 1:Mdl_T_amb_test.P;            % Contains 2 observations
% estsample_T_amb_test = (Mdl_T_amb_test.P + 1):Size_T_amb_test; % Contains the remaining observations
% EstMdl_H_Env = estimate(Mdl_T_amb_test, Raw.T_amb_test(estsample_T_amb_test), 'Y0', Raw.T_amb_test(presample_T_amb_test));
% p.a_T_amb_test    = cell2mat(EstMdl_T_amb_test.AR(1)); % AR coefficient, alpha
% p.w_T_amb_test    = EstMdl_T_amb_test.Variance;        % Variance which is rooted to obtain error
% p.c_T_amb_test    = EstMdl_T_amb_test.Constant;        % Constant or intercept
% g.T_amb_test(1) = Raw.T_amb_test(1);
% for i = 1:51086           
%     e.T_amb_test(i)   = randn*sqrt(p.w_H_Env);
%     g.T_amb_test(i+1) = p.c_T_amb_test + p.a_T_amb_test*g.T_amb_test(i) + e.T_amb_test(i);    % Generated data used to test the model in simulation
% end
% u.T_amb_generated_test = griddedInterpolant(t, g.T_amb_test');
% days = 86400;
% plot(t/days, Raw.T_amb_test, 'r', t/days, u.T_amb_generated_test(t)', 'b')
% legend('Raw Data', 'AR Model')
% xlabel('Time (days)')
% ylabel('T_A_n_t_i_S_u_r_g_e (^oC)')
% ax = gca;
% ax.FontSize = p.font_size;


% ARIMA MODEL FOR v.n
Raw.v  = sum(u.s(t),2);
Size_s = size(Raw.v,1);
t1 = t(1:51086,1);
%Mdl_s = arima(1,0,0); % First-order Auto-Regressive Model plus constant
dep = 0.998; % Dependence of the next L measurement on the previous one
Mdl_s = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant
presample_s = 1:Mdl_s.P;               % Contains 2 observations
estsample_s = (Mdl_s.P + 1):Size_s; % Contains the remaining observations
EstMdl_s = estimate(Mdl_s, Raw.v(estsample_s), 'Y0', Raw.v(presample_s));
p.a_s = cell2mat(EstMdl_s.AR(1)); % AR coefficient, alpha
p.w_s = EstMdl_s.Variance;        % Variance which is rooted to obtain error
p.c_s = EstMdl_s.Constant;        % Constant or intercept
g.s(1)  = Raw.v(1);  % Specify the inital value of the inlet flowrate
for i = 1:51086           
    e.s(i)    = randn*sqrt(p.w_s);
    g.s(i+1)  = p.c_s + p.a_s*g.s(i) + e.s(i);    % Generated data used to test the model in simulation
end


% % % ARIMA MODEL FOR u_s
Raw.s = u.s(t);
% Size_s = size(Raw.s,1);
% t1 = t(1:51086,1);
% 
% % Create Model Template
% %Mdl_s = arima(1,0,0); 
% dep = 0.999; % Dependence of the next L measurement on the previous one
% Mdl_s = arima('AR',{dep}); % First-order Auto-Regressive Model plus constant
% 
% presample_s = 1:Mdl_s.P;            % Contains 2 observations
% estsample_s = (Mdl_s.P + 1):Size_s; % Contains the remaining observations
% 
% for i = 1:5
%     EstMdl_s(i) = estimate(Mdl_s, Raw.s(estsample_s,i), 'Y0', Raw.s(presample_s,i));
%     p.a_s(i) = cell2mat(EstMdl_s(i).AR(1)); % AR coefficient, alpha
%     p.w_s(i) = EstMdl_s(i).Variance;        % Variance which is rooted to obtain error
%     p.c_s(i) = EstMdl_s(i).Constant;        % Constant or intercept
%     % AR Model Predictions
%     k.s(1,i)  = Raw.s(1,i);  % Specify the inital value of the inlet flowrate
% end
% 
% for j = 1:51086 % From here onwards, make predictions
%     e.s(j,:) = randn.*sqrt(p.w_s);
%     k.s(j+1,:) = p.c_s + p.a_s.*k.s(j,:) + e.s(j,:);
% end

g_s = g.s';
yval_g = [];
Size_g_s = size(g_s,1);
ipt_g = findchangepts(g_s,'Statistic','std','MinThreshold',100);
intvls_g = [1 ipt_g' length(g_s)];
for k1_g = 1:length(intvls_g)-2
    interval_means_g(k1_g) = mean(g_s(intvls_g(k1_g):intvls_g(k1_g+1)-1));
end
interval_means_g(k1_g+1) = mean(g_s(intvls_g(k1_g):intvls_g(end)));
xval_g = [];
yval_g = [];
for k1_g = 1:length(intvls_g)-1
    xval_g = intvls_g(k1_g):intvls_g(k1_g+1);
ynext_g = interval_means_g(k1_g)*ones(1,(length(xval_g)-1));
yval_g = [yval_g ynext_g];
end
yval_g(end+1) = 0;
g.s = round(yval_g');
g.s = g.s + 1;
g.s = min(5,g.s);


k_s = zeros(length(t),5);
for i = 1:length(g.s)
    for j = 1:g.s(i)
        k_s(i,j) = 1;
    end
end
u.s_generated = griddedInterpolant(t, k_s, "previous"); 


% for i = 1:5
%     figure(i)
%     plot(t/86400, Raw.s(:,i), 'r', t/86400, u.s_generated.Values(:,i)', 'b')
%     legend('Raw Data', 'AR Model')
%     xlabel('Time (days)')
%     ylabel('s')
%     ax = gca;
%     ax.FontSize = p.font_size;
%     ylim([-0.5 1.5]);
% end
% 

    % plot(t/86400, g.s, t/86400, Raw.v)
    % legend('Raw Data', 'AR Model')
    % xlabel('Time (days)')
    % ylabel('Number of RPs On')
    % ax = gca;
    % ax.FontSize = p.font_size;
    % ax.YAxis(1).Exponent = 0;


end