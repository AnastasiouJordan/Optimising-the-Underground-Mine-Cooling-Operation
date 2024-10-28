%% System of ODEs: Pre-cooling Towers
%  Jordan Anastasiou, 2022-06
%  This code is for the pre-cooling towers with both towers
%  being grouped together to form one system over which the
%  mass and energy balances were performed. There is a flowrate
%  of water as well as a flowrate of air through the tower.
clc
clear
clf


%% Define parameters
p.regressedparameterfields = {'m_Air'}; 

p.rho_Water = 1000;      % kg/m3,  Density of water
p.rho_Air   = 1.225;     % kg/m3,  Density of air
p.H_outAir  = 100;       % %,      Humidity of air leaving towers
p.h_0       = 0.10186;   % kJ/kg,  Reference specific enthalpy
p.T_0       = 0.01;      % oC,     Reference temperature
p.A         = 1.01;      % -,      Enthalpy equation constant
p.B         = 1.89;      % -,      Enthalpy equation constant
p.C         = 2501;      % -,      Enthalpy equation constant
p.C1        = -5.8E+3;   % -,      Pressure equation constant
p.C2        = 1.391;     % -,      Pressure equation constant
p.C3        = -4.864E-2; % -,      Pressure equation constant
p.C4        = 4.176E-5;  % -,      Pressure equation constant
p.C5        = -1.445E-8; % -,      Pressure equation constant
p.C6        = 6.546;     % -,      Pressure equation constant
p.M_r       = 0.62198;   % -,      Ratio between molar mass water and air
p.P_Tot     = 101.325;   % kPa,    Atmospheric pressure
p.C_p       = 4.1831;    % kJ/kgC, Heat capacity of water

p.height_PT = 3;         % m, Height of the Pre-cooling Tower basins
p.length_PT = 14.8;      % m, Length of the Pre-cooling Tower basins
p.width_PT  = 10.45*2;   % m, Width of the Pre-cooling Tower basins
p.area_PT   = p.length_PT*p.width_PT; % m2, Area of the Pre-cooling Tower basins
p.volume_PT = p.area_PT*p.height_PT;  % m3, Volume of the Pre-cooling Tower basins

p.m_PTAmax  = p.volume_PT/2*1000;     % kg, Maximum mass capacity of PT A basin
p.m_PTBmax  = p.volume_PT/2*1000;     % kg, Maximum mass capacity of PT B basin

% Initial guess for m_Air
p.m_Air     = 700;       % kg/s, Estimated mass flowrate of air.
                         % This is based on a L/G ratio of 0.5 from 
                         % literature with the water flowrate through
                         % the towers having an average range of 
                         % 300 - 500 kg/s.
pmAirVec = S2V(p, p.regressedparameterfields); % Convert the unknown parameter to a vector

%% Define exogeneous inputs

load SavedInterpolantsPT.mat


%% Define state structure and initial conditions
s.statefields = {'m_PT','h_PT'}; % Field names for each state

x0.m_PT = u.L_PTA(0).*(p.m_PTAmax + p.m_PTBmax)/100;
x0.h_PT = 80;                     % kJ/kgK, Initial value for enthalpy
                                  % of water leaving PTs
x0_vec  = S2V(x0, s.statefields); % Convert initial state values to vector

%% Load Kalman Filtered Values

load KalmanFilterPT.mat

%% Simulate system of ODEs
[~, x_vec] = ode45(@(t, x) PreCoolingTowersODEs(s, p, x, u, t), t, x0_vec);
x = V2S(x_vec', s.statefields);
v = PTIntermediates(x, u, p, t);
save PreCoolingTowers.mat v u

%% Plot
% State variables. Plot dm_PT/dt and dh_PT/dt calculated using ODE.
font_size = 18;

figure (1)
title('Prediction Results with Assumed Parameter Values');
subplot(2,1,1)
plot(t/86400, u.T_PT(t), t/86400, v.T_PT, 'LineWidth', 1.5);  
legend('measured', 'predicted', 'FontSize', font_size);
%text(3100000, 2, sprintf('mAir = %.2f kJ/Ks', p.m_Air));
xlabel('Time (days)');
ylabel('T_P_T (^oC)');
ax = gca;
ax.FontSize = font_size;
subplot(2,1,2)
plot(t/86400, u.L_PTA(t), t/86400, v.L_PT, 'LineWidth', 1.5);
legend('measured', 'predicted', 'FontSize', font_size);
%text(3100000, 2, sprintf('mAir = %.2f kJ/Ks', p.m_Air));
xlabel('Time (days)');
ylabel('L_P_T (%)');
ax = gca;
ax.FontSize = font_size;


%% Regression

% options = optimoptions('lsqnonlin', 'StepTolerance', 1e-9,...
%                        'Algorithm','trust-region-reflective',...
%                        'FiniteDifferenceType','central',...
%                        'TypicalX', 650,...
%                        'Display', 'iter-detailed');
% 
% p_est    = lsqnonlin(@(pmAirVec) PTCalcError(pmAirVec, u, p, s, t), pmAirVec, 100, 900, options)
% 
% [E, x, v] = PTCalcError(p_est, u, p, s, t);
% 
% %Plot results using parameter estimate/regressed parameter
% figure (2)
% title('Prediction Results with Regressed Parameter Values');
% subplot(2,1,1)
% plot(t, u.T_PT(t), t, v.T_PT);  
% legend('measured', 'predicted');
% text(3100000, 2, sprintf('mAir = %.2f kJ/Ks', p.m_Air));
% xlabel('Time (s)');
% ylabel('T_P_T (^oC)');
% subplot(2,1,2)
% plot(t, u.L_PTA(t), t, v.L_PT);
% legend('measured', 'predicted');
% text(3100000, 2, sprintf('mAir = %.2f kJ/Ks', p.m_Air));
% xlabel('Time (s)');
% ylabel('L_P_T (%)');
% ax = gca;
% ax.FontSize = font_size;

%% Likelihood Profiles

% SSR_Level = []; % Sum of squared residuals for level
% SSR_Temp  = []; % Sum of squared residuals for temperature
% soln_space = 0.01:200:3000;
% for pmAirvec = soln_space
%     [E, x, v] = PTCalcError(pmAirvec, u, p, s, t);
%     E_2 = E.^2;
%     sum_E_2_Level = sum(E_2(:,1));
%     sum_E_2_Temp  = sum(E_2(:,2));
%     SSR_Level = [SSR_Level sum_E_2_Level];
%     SSR_Temp  = [SSR_Temp sum_E_2_Temp];
% end
% 
% figure(3)
% title('Likelihood Profile for parameter mAir and Level');
% LLRatio_Level = 2*log(SSR_Level/min(SSR_Level));
% plot(soln_space, LLRatio_Level);
% hold on
% yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
% hold off
% xlabel('m_A_i_r (m/s)');
% ylabel('Negative Log Likelihood Ratio');
% xlim([0 soln_space(end)]);
% x1 = interp1(LLRatio_Level, soln_space, 0);
% zero_point = find(soln_space == x1);
% x2 = interp1(LLRatio_Level(1:zero_point), soln_space(1:zero_point), 2.71);
% x3 = interp1(LLRatio_Level(zero_point:end), soln_space(zero_point:end), 2.71);
% xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
% xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
% xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
% ax = gca;
% ax.FontSize = font_size;
% 
% figure(4)
% LLRatio_Temp = 2*log(SSR_Temp/min(SSR_Temp));
% plot(soln_space, LLRatio_Temp);
% hold on
% yline(2.71,'-',{'Chi-Square Threshold'});
% hold off
% xlabel('m_A_i_r  (m/s)');
% ylabel('Negative Log Likelihood Ratio');
% xlim([0 3000]);
% x1 = interp1(LLRatio_Temp, soln_space, 0);
% zero_point = find(soln_space == x1);
% x2 = interp1(LLRatio_Temp(1:zero_point), soln_space(1:zero_point), 2.71);
% x3 = interp1(LLRatio_Temp(zero_point:end), soln_space(zero_point:end), 2.71);
% xline(x1, '--', {'Optimal Parameter Value'});
% xline(x2, '--', {'90% Confidence Interval'});
% xline(x3, '--', {'90% Confidence Interval'});
% ax = gca;
% ax.FontSize = font_size;
% 
% 
% %% MAPE
% forecast1 = v.T_PT;
% observed1 = u.T_PT(t);
% ABS1 = abs((observed1 - forecast1)./forecast1)*100;
% MAPE1 = 1/length(t)*sum(ABS1);
% 
% forecast2 = v.L_PT;
% observed2 = u.L_PTA(t);
% ABS2 = abs((forecast2 - observed2)./observed2)*100;
% MAPE2 = 1/length(t)*sum(ABS2);
% 
% MAPEPT = [MAPE1 MAPE2]