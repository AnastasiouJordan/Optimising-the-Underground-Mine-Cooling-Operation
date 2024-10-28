%% System of ODEs: Chill Dam
%  Jordan Anastasiou, 07-2022
%  This code is for the chill dam

clc
clear
clf

%% Define parameters
p.regressedparameterfields = {'UA_CD'};

p.rho_Water = 1000;    % kg/m3,  Density of water
p.h_0       = 0.10186; % kJ/kg,  Reference specific enthalpy
p.T_0       = 0.01;    % oC,     Reference temperature
p.C_p       = 4.1831;  % kJ/kgK, Heat capacity of water
p.m_CDmax   = 8000000; % kg,     Maximum mass capacity of CD


p.height_CD = 6.7; % m, Height of the Chill Dam
p.length_CD = 48;  % m, Length of the Chill Dam
p.width_CD  = 25;  % m, Width of the Chill Dam

p.area_CD   = p.length_CD*p.width_CD; % m2, Area of the Chill Dam

% Initial guess for UA_CD

p.UA_CD     = 190; % kJ/Ks,  Heat transfer coefficient for the rate
                   % of heat transfer between the water in the CD and
                  % the environment through the dam walls

pUAvec = S2V(p, p.regressedparameterfields); % Convert the unknown parameter to a vector

%% Define exogeneous inputs

load SavedInterpolantsCD.mat

%% Define state structure and initial conditions

s.statefields = {'m_CD','h_CD'};   % Field names for each state  
x0.m_CD = u.L_CD(0)*p.m_CDmax/100;
x0.h_CD = 5;                      % kJ/kgK, Initial value for enthalpy
                                   % of water leaving CD
x0_vec = S2V(x0, s.statefields);

%% Load Kalman Filtered Values

load KalmanFilterCD.mat

%% Simulate system of ODEs

[~, x_vec] = ode45(@(t, x) ChillDamODEs(s, p, x, u, t), t, x0_vec);
x = V2S(x_vec', s.statefields);
v = CDIntermediates(x, u, p, t);

%% Plot
% Plot results using initial assumption for the parameter value
font_size = 18;
% figure (1)
% title('Prediction Results with Assumed Parameter Values');
% subplot(2,1,1)
% plot(t/86400, u.T_CD(t), t/86400, v.T_CD);  
% legend('measured', 'predicted','FontSize',font_size);
% %text(3100000, 2, sprintf('UA = %.2f kJ/Ks', p.UA_CD));
% xlabel('Time (days)');
% ylabel('T_C_D (^oC)');
% ax = gca;
% ax.FontSize = font_size;
% subplot(2,1,2)
% plot(t/86400, u.L_CD(t), t/86400, v.L_CD);
% legend('measured', 'predicted','FontSize',font_size);
% xlabel('Time (days)');
% ylabel('L_C_D (%)');
% ax = gca;
% ax.FontSize = font_size;

figure(3)
font_size = 18;
plot(t/86400, u.T_CD(t), t/86400, x.h_CD');  
legend('measured', 'predicted','FontSize',font_size);
%text(3100000, 2, sprintf('UA = %.2f kJ/Ks', p.UA_CD));
xlabel('Time (days)');
ylabel('T_C_D (^oC)');
ax = gca;
ax.FontSize = font_size;

figure(4)
plot(t/86400, u.L_CD(t), t/86400, v.L_CD);
legend('measured', 'predicted','FontSize',font_size);
xlabel('Time (days)');
ylabel('L_C_D (%)');
ax = gca;
ax.FontSize = font_size;


%% Regression

options = optimoptions('lsqnonlin', 'StepTolerance', 1e-9,...
                       'Algorithm','trust-region-reflective',...
                       'FiniteDifferenceType','central',...
                       'TypicalX', 200,...
                       'Display', 'iter-detailed');

p_est    = lsqnonlin(@(pUAvec) CDCalcError(pUAvec, u, p, s, t), pUAvec, 50, 250, options)


[E, x, v] = CDCalcError(p_est, u, p, s, t);


% Plot results using parameter estimate/regressed parameter

figure (2)
title('Prediction Results with Regressed Parameter Values');
subplot(2,1,1)
plot(t, u.T_CD(t), t, v.T_CD);  
legend('measured', 'predicted');
text(3100000, 2, sprintf('UA = %.2f kJ/Ks', p.UA_CD));
xlabel('Time (s)');
ylabel('T_C_D (^oC)');
subplot(2,1,2)
plot(t, u.L_CD(t), t, v.L_CD);
legend('measured', 'predicted');
text(3100000, 2, sprintf('UA = %.2f kJ/Ks', p.UA_CD));
xlabel('Time (s)');
ylabel('L_C_D (%)');
ax = gca;
ax.FontSize = font_size;

%% Likelihood profile
SSR_Level = []; % Sum of squared residuals for level
SSR_Temp  = []; % Sum of squared residuals for temperature
soln_space = 0.01:10:1000;
for pUAvec = soln_space
    [E, x, v] = CDCalcError(pUAvec, u, p, s, t);
    E_2 = E.^2;
    sum_E_2_Level = sum(E_2(:,1));
    sum_E_2_Temp  = sum(E_2(:,2));
    SSR_Level = [SSR_Level sum_E_2_Level];
    SSR_Temp  = [SSR_Temp sum_E_2_Temp];
end

figure(3)
title('Likelihood Profile for parameter UA and Level');
LLRatio_Level = 2*log(SSR_Level/min(SSR_Level));
plot(soln_space, LLRatio_Level);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
hold off
xlabel('UA_C_D (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([0 1000]);
x1 = interp1(LLRatio_Level, soln_space, 0);
zero_point = find(soln_space == x1);
x2 = interp1(LLRatio_Level(1:zero_point), soln_space(1:zero_point), 2.71);
x3 = interp1(LLRatio_Level(zero_point:end), soln_space(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size;

figure(4)
LLRatio_Temp = 2*log(SSR_Temp/min(SSR_Temp));
plot(soln_space, LLRatio_Temp);
hold on
yline(2.71,'-',{'Chi-Square Threshold'});
hold off
xlabel('UA_C_D (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([0 1000]);
x1 = interp1(LLRatio_Temp, soln_space, 0);
zero_point = find(soln_space == x1);
x2 = interp1(LLRatio_Temp(1:zero_point), soln_space(1:zero_point), 2.71);
x3 = interp1(LLRatio_Temp(zero_point:end), soln_space(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'});
xline(x2, '--', {'90% Confidence Interval'});
xline(x3, '--', {'90% Confidence Interval'});
ax = gca;
ax.FontSize = font_size;

%% MAPE for all units
% CD
forecastTCD = v.T_CD;
observedTCD = u.T_CD(t);
ABSTCD = abs((forecastTCD - observedTCD)./observedTCD)*100;
MAPETCD = 1/length(t)*sum(ABSTCD);

forecastLCD = v.L_CD;
observedLCD = u.L_CD(t);
ABSLCD = abs((forecastLCD - observedLCD)./observedLCD)*100;
MAPELCD = 1/length(t)*sum(ABSLCD);

MAPECD = [MAPETCD MAPELCD];

%% PT
load PreCoolingTowers.mat
forecastTPT = v.T_PT;
observedTPT = u.T_PT(t);
ABSTPT = abs((observedTPT - forecastTPT)./forecastTPT)*100;
MAPETPT = 1/length(t)*sum(ABSTPT);

forecastLPT = v.L_PT;
observedLPT = u.L_PTA(t);
ABSLPT = abs((forecastLPT - observedLPT)./observedLPT)*100;
MAPELPT = 1/length(t)*sum(ABSLPT);

MAPEPT = [MAPETPT MAPELPT];

%% RPs
load FridgePlants.mat
N = 400:650;
% RP1
forecast1 = v.T_RP(:,1);
observed1 = u.T_RP.Values(N,1);
ABS1 = abs((forecast1 - observed1)./observed1)*100;
MAPE1 = 1/length(N)*sum(ABS1);

% RP2
forecast2 = v.T_RP(:,2);
observed2 = u.T_RP.Values(N,2);
ABS2 = abs((forecast2 - observed2)./observed2)*100;
MAPE2 = 1/length(N)*sum(ABS2);

% RP3
forecast3 = v.T_RP(:,3);
observed3 = u.T_RP.Values(N,3);
ABS3 = abs((forecast3 - observed3)./observed3)*100;
MAPE3 = 1/length(N)*sum(ABS3);

% RP4
forecast4 = v.T_RP(:,4);
observed4 = u.T_RP.Values(N,4);
ABS4 = abs((forecast4 - observed4)./observed4)*100;
MAPE4 = 1/length(N)*sum(ABS4);

% RP5
forecast5 = v.T_RP(:,5);
observed5 = u.T_RP.Values(N,5);
ABS5 = abs((forecast5 - observed5)./observed5)*100;
MAPE5 = 1/length(N)*sum(ABS5);

MAPERP = [MAPE1 MAPE2 MAPE3 MAPE4 MAPE5];

%% SD
load SquareDam.mat
forecastLSD = v.L_SD;
observedLSD = u.L_SD(t);
ABSLSD = abs((observedLSD - forecastLSD)./forecastLSD)*100;
MAPELSD = 1/length(t)*sum(ABSLSD);


%% Plot
ALL_ABS = [ABSLSD; ABSLPT; ABSTPT; ABS1; ABS2; ABS3; ABS4; ABS5; ABSLCD; ABSTCD];
group = [repmat({'SD Level'},51087,1); repmat({'PT Level'},51087,1); repmat({'PT Temp'},51087,1); repmat({'RP1 Temp'},251,1);...
    repmat({'RP2 Temp'},251,1);repmat({'RP3 Temp'},251,1);repmat({'RP4 Temp'},251,1);repmat({'RP5 Temp'},251,1);repmat({'CD Level'},51087,1);repmat({'CD Temp'},51087,1)];
bh=boxplot(ALL_ABS,group,'PlotStyle','compact','MedianStyle','line','BoxStyle','outline','Symbol','.r','LabelOrientation','horizontal');
% set(gca,'TickLabelInterpreter','tex');
% set(gca, 'xticklabel', {'L_S_D', 'L_P_T', 'T_P_T', 'T_R_P_1', 'T_R_P_2','T_R_P_3','T_R_P_4', 'T_R_P_5', 'L_C_D', 'T_C_D'})
ylim([-2 85]);
ylabel('Absolute Percentage Error');
ax = gca;
ax.FontSize = 30;
set(bh,'LineWidth', 2);
