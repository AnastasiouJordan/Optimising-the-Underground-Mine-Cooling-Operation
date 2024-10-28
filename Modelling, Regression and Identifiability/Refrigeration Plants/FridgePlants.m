%% System of ODEs: Refrigeration Plants
%  Jordan Anastasiou, 2022-03-15
%  This code is for the refrigeration plants where
%  each plant is modelled as a system, j.

clc
clear
clf


%% Define parameters
p.rho_Water = 1000;      % kg/m3,  Density of water
p.m_RPj     = 10000;     % kg,     Mass held by each fridge plant
                         %         Differs from v.m_RPj which is the mass
                         %         flowrate through each plant respectively
p.h_0       = 0.10186;   % kJ/kg,  Reference specific enthalpy
p.T_0       = 0.01;      % oC,     Reference temperature
p.C_p       = 4.1831;    % kJ/kgC, Heat capacity of water

% Initial guesses for unknown parameters UA_RPj and UA_amb
p.UA_RP1  = 10000;       % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 1
p.UA_RP2  = 15000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 2
p.UA_RP3  = 5000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 3
p.UA_RP4  = 32500;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 4
p.UA_RP5  = 58000;        % kJ/Ks,  Evaporator heat transfer coefficient for fridge plant 5
p.UA_RP   = [p.UA_RP1, p.UA_RP2, p.UA_RP3, p.UA_RP4, p.UA_RP5];

p.UA_amb1 = -10;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 1
p.UA_amb2 = -10;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 2
p.UA_amb3 = -30;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 3
p.UA_amb4 = -30;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 4
p.UA_amb5 = -50;         % kJ/Ks,  Ambient heat transfer coefficient for fridge plant 5
p.UA_amb  = [p.UA_amb1, p.UA_amb2, p.UA_amb3, p.UA_amb4, p.UA_amb5];
p.regressedparameterfields = {'UA_RP', 'UA_amb'}; % Parameter structure field names

pUAvec = S2V(p, p.regressedparameterfields); % Convert the unknown parameters to a vector of unknown parameters

%% Define exogeneous inputs

load SavedInterpolantsRP.mat
load ArimaModelsRP.mat

%% Simulate system of ODEs
%N = length(t); % Simulate only a piece of the data. N = length(t) to simulate full set.
N = 51087;
h_RP0 = p.C_p * (u.T_RP(0) - p.T_0) + p.h_0;
options = odeset('OutputFcn', @odeplot);
tic
[~, h_RP] = ode45(@(t, h_RP) FridgePlantsODEs(p, h_RP, u, t), t(1:N), h_RP0, options);
toc
v = RPIntermediates(h_RP', u, p, t(1:N));

save FridgePlants.mat v u

%% Plot

% Measured vs calculated. Plot the variables for which there is both
% measured and calculated data. These are also the variables used in
% the objective function.
figure (1)
font_size = 18;
plot(t(1:N)/86400, u.T_RP.Values(1:N,1), t(1:N)/86400, v.T_RP(:,1),t(1:N)/86400, u.s.Values(1:N,1), 'm-'); 
legend('measured', 'predicted','on/off status', 'FontSize', font_size);
xlabel('Time (days)');
ylabel('T_R_P_1 (^oC)');
ax = gca;
ax.FontSize = font_size;
figure(2)
plot(t(1:N)/86400, u.T_RP.Values(1:N,2), t(1:N)/86400, v.T_RP(:,2),t(1:N)/86400, u.s.Values(1:N,2), 'm-');
legend('measured', 'predicted','on/off status', 'FontSize', font_size);
xlabel('Time (days)');
ylabel('T_R_P_2 (^oC)');
ax = gca;
ax.FontSize = font_size;
figure(3)
plot(t(1:N)/86400, u.T_RP.Values(1:N,3), t(1:N)/86400, v.T_RP(:,3),t(1:N)/86400, u.s.Values(1:N,3), 'm-');
legend('measured', 'predicted','on/off status', 'FontSize', font_size);
xlabel('Time (days)');
ylabel('T_R_P_3 (^oC)');
ax = gca;
ax.FontSize = font_size;
figure(4)
plot(t(1:N)/86400, u.T_RP.Values(1:N,4), t(1:N)/86400, v.T_RP(:,4),t(1:N)/86400, u.s.Values(1:N,4), 'm-');
legend('measured', 'predicted','on/off status', 'FontSize', font_size);
xlabel('Time (days)');
ylabel('T_R_P_4 (^oC)');
ax = gca;
ax.FontSize = font_size;
figure(5)
plot(t(1:N)/86400, u.T_RP.Values(1:N,5), t(1:N)/86400, v.T_RP(:,5),t(1:N)/86400, u.s.Values(1:N,5), 'm-');
legend('measured', 'predicted','on/off status', 'FontSize', font_size);
xlabel('Time (days)');
ylabel('T_R_P_5 (^oC)');
ax = gca;
ax.FontSize = font_size;


figure (6)
subplot(5,1,1)
plot(t(1:N), u.F_inRP.Values(1:N,1), t(1:N), v.m_RP(1:N,1),t(1:N), u.s.Values(1:N,1)*200); 
legend('measured', 'predicted','on/off status');
xlabel('Time (s)');
ylabel('F_R_P_1 (L/s)');
subplot(5,1,2)
plot(t(1:N), u.F_inRP.Values(1:N,2), t(1:N), v.m_RP(:,2),t(1:N), u.s.Values(1:N,2)*200); 
legend('measured', 'predicted','on/off status');
xlabel('Time (s)');
ylabel('F_R_P_2 (L/s)');
subplot(5,1,3)
plot(t(1:N), u.F_inRP.Values(1:N,3), t(1:N), v.m_RP(:,3),t(1:N), u.s.Values(1:N,3)*200); 
legend('measured', 'predicted','on/off status');
xlabel('Time (s)');
ylabel('F_R_P_3 (L/s)');
subplot(5,1,4)
plot(t(1:N), u.F_inRP.Values(1:N,4), t(1:N), v.m_RP(:,4),t(1:N), u.s.Values(1:N,4)*200); 
legend('measured', 'predicted','on/off status');
xlabel('Time (s)');
ylabel('F_R_P_4 (L/s)');
subplot(5,1,5)
plot(t(1:N), u.F_inRP.Values(1:N,5), t(1:N), v.m_RP(:,5),t(1:N), u.s.Values(1:N,5)*200); 
legend('measured', 'predicted','on/off status');
xlabel('Time (s)');
ylabel('F_R_P_5 (L/s)');


%% Regression

% options = optimoptions('lsqnonlin', 'StepTolerance', 1e-7,...
%                        'Algorithm','trust-region-reflective',...
%                        'FiniteDifferenceType','central',...
%                        'Display', 'iter-detailed');
% 
% p_est1    = lsqnonlin(@(pUAvec) RPCalcError1(pUAvec, u, p, t, N), pUAvec(:,1), [1000; 100], [10000; 500], options);
% p_est2    = lsqnonlin(@(pUAvec) RPCalcError2(pUAvec, u, p, t, N), pUAvec(:,2), [1000; 100], [10000; 500], options);
% p_est3    = lsqnonlin(@(pUAvec) RPCalcError3(pUAvec, u, p, t, N), pUAvec(:,3), [1000; 100], [10000; 500], options);
% p_est4    = lsqnonlin(@(pUAvec) RPCalcError4(pUAvec, u, p, t, N), pUAvec(:,4), [1000; 100], [25000; 500], options);
% p_est5    = lsqnonlin(@(pUAvec) RPCalcError5(pUAvec, u, p, t, N), pUAvec(:,5), [1000; 100], [25000; 500], options);
% 
% p_est = [p_est1; p_est2; p_est3; p_est4; p_est5];
% 
% 
% [E1, h_RP, v] = RPCalcError1(p_est1, u, p, t, N);
% [E2, h_RP, v] = RPCalcError2(p_est2, u, p, t, N);
% [E3, h_RP, v] = RPCalcError3(p_est3, u, p, t, N);
% [E4, h_RP, v] = RPCalcError4(p_est4, u, p, t, N);
[E5, h_RP, v] = RPCalcError5(p_est5, u, p, t, N);


%Plot results using parameter estimate/regressed parameter
% figure (3)
% font_size = 17;
% subplot(5,1,1)
% plot(t(1:N), u.T_RP.Values(1:N,1), t(1:N), v.T_RP(:,1),t(1:N), u.s.Values(1:N,1), 'm-'); 
% legend('measured', 'predicted','on/off status', 'FontSize', font_size);
% xlabel('Time (s)');
% ylabel('T_R_P_1 (^oC)');
% ax = gca;
% ax.FontSize = font_size;
% subplot(5,1,2)
% plot(t(1:N), u.T_RP.Values(1:N,2), t(1:N), v.T_RP(:,2),t(1:N), u.s.Values(1:N,2), 'm-');
% legend('measured', 'predicted','on/off status', 'FontSize', font_size);
% xlabel('Time (s)');
% ylabel('T_R_P_2 (^oC)');
% ax = gca;
% ax.FontSize = font_size;
% subplot(5,1,3)
% plot(t(1:N), u.T_RP.Values(1:N,3), t(1:N), v.T_RP(:,3),t(1:N), u.s.Values(1:N,3), 'm-');
% legend('measured', 'predicted','on/off status', 'FontSize', font_size);
% xlabel('Time (s)');
% ylabel('T_R_P_3 (^oC)');
% ax = gca;
% ax.FontSize = font_size;
% subplot(5,1,4)
% plot(t(1:N), u.T_RP.Values(1:N,4), t(1:N), v.T_RP(:,4),t(1:N), u.s.Values(1:N,4), 'm-');
% legend('measured', 'predicted','on/off status', 'FontSize', font_size);
% xlabel('Time (s)');
% ylabel('T_R_P_4 (^oC)');
% ax = gca;
% ax.FontSize = font_size;
% subplot(5,1,5)
% plot(t(1:N), u.T_RP.Values(1:N,5), t(1:N), v.T_RP(:,5),t(1:N), u.s.Values(1:N,5), 'm-');
% legend('measured', 'predicted','on/off status', 'FontSize', font_size);
% xlabel('Time (s)');
% ylabel('T_R_P_5 (^oC)');
% ax = gca;
% ax.FontSize = font_size;


%% Likelihood Profiles

% The following code takes the vector of parameters, pUAvec, which is
% a 2x5 vector structured as follows:
% pUA_RP1  pUA_RP2  pUA_RP3  pUA_RP4  pUA_RP5
% pUA_amb1 pUA_amb2 pUA_amb3 pUA_amb4 pUA_amb5
% and it varies only ONE parameter for each fridge plant at a time.
% In other words, for fridge plant 1, pUA_RP1 will be varied over a
% specified solution space, while pUA_amb1 will be kept constant, 
% and the parameters of fridge plant 2 through 5 are also kept constant.
% The purpose of this is to determine the effect of varying each parameter
% on the error in the temperature prediction.
% The parameter being varied is indicated by a position in the overall
% vector of parameters, by specifying 'varyRowIndex' and 'varyColumnIndex'.

font_size = 30; % Define the font size to be used in all plotted figures

% FRIDGE PLANT 1
% For pUA_RP1
varyRowIndex    = 1; % Define the row position of the parameter you want to change
varyColumnIndex = 1; % Define the column position of the parameter you want to change

SSR_11 = [];                    % Pre-define a sum of squared residuals
soln_space11 = 8000:1000:23000; % Define the solution space for parameter at position 1,1 in pUAvec

for soln_space = soln_space11     % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E1, x, v] = RPCalcError1(modifiedVector, u, p, t, N); % Calculate the error at each value
                E11_2 = E1.^2;               % Square the error at each value
                sum_E11_2 = sum(E11_2);      % Sum the squared error at each value
                SSR_11 = [SSR_11 sum_E11_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end

figure(4)
subplot(2,1,1)
title('Likelihood Profile for RP1, parameter UA_RP1');
LLRatio_UARP1 = 2*log(SSR_11/min(SSR_11));
plot(soln_space11, LLRatio_UARP1);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
hold off
xlabel('UA_R_P_1 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([soln_space11(1) soln_space11(end)]);
x1 = interp1(LLRatio_UARP1, soln_space11, 0);
zero_point = find(soln_space11 == x1);
x2 = interp1(LLRatio_UARP1(1:zero_point), soln_space11(1:zero_point), 2.71);
x3 = interp1(LLRatio_UARP1(zero_point:end), soln_space11(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

% For pUA_amb1
varyRowIndex    = 2; % Define the row position of the parameter you want to change
varyColumnIndex = 1; % Define the column position of the parameter you want to change

SSR_21 = [];                 % Pre-define a sum of squared residuals
soln_space21 = 0:100:400; % Define the solution space for parameter at position 2,1 in pUAvec

for soln_space = soln_space21     % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E1, x, v] = RPCalcError1(modifiedVector, u, p, t, N); % Calculate the error at each value
                E21_2 = E1.^2;               % Square the error at each value
                sum_E21_2 = sum(E21_2);      % Sum the squared error at each value
                SSR_21 = [SSR_21 sum_E21_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end

subplot(2,1,2)
title('Likelihood Profile for RP1, parameter UA_amb1');
LLRatio_UAamb1 = 2*log(SSR_21/min(SSR_21));
plot(soln_space21, LLRatio_UAamb1);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
hold off
xlabel('UA_a_m_b_1 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([soln_space21(1) soln_space21(end)]);
x1 = interp1(LLRatio_UAamb1, soln_space21, 0);
zero_point = find(soln_space21 == x1);
x2 = interp1(LLRatio_UAamb1(1:zero_point), soln_space21(1:zero_point), 2.71);
x3 = interp1(LLRatio_UAamb1(zero_point:end), soln_space21(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

% FRIDGE PLANT 2
% For pUA_RP2
%N = 38096:38098;
varyRowIndex    = 1; % Define the row position of the parameter you want to change
varyColumnIndex = 2; % Define the column position of the parameter you want to change

SSR_12 = [];                    % Pre-define a sum of squared residuals
soln_space12 = -10000:50000:500000; % Define the solution space for parameter at position 1,2 in pUAvec

for soln_space = soln_space12     % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E2, x, v] = RPCalcError2(modifiedVector, u, p, t, N); % Calculate the error at each value
                E12_2 = E2.^2;               % Square the error at each value
                sum_E12_2 = sum(E12_2);      % Sum the squared error at each value
                SSR_12 = [SSR_12 sum_E12_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end

figure(5)
subplot(2,1,1)
title('Likelihood Profile for RP2, parameter UA_RP2');
LLRatio_UARP2 = 2*log(SSR_12/min(SSR_12));
plot(soln_space12, LLRatio_UARP2);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
ax = gca;
ax.FontSize = font_size; 
hold off
xlabel('UA_R_P_2 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([soln_space12(1) soln_space12(end)]);
x1 = interp1(LLRatio_UARP2, soln_space12, 0);
zero_point = find(soln_space12 == x1);
x2 = interp2(LLRatio_UARP2(1:zero_point), soln_space12(1:zero_point), 2.71);
x3 = interp2(LLRatio_UARP2(zero_point:end), soln_space12(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

% For pUA_amb2
varyRowIndex    = 2; % Define the row position of the parameter you want to change
varyColumnIndex = 2; % Define the column position of the parameter you want to change

SSR_22 = [];                 % Pre-define a sum of squared residuals
soln_space22 = -5000:100:10000; % Define the solution space for parameter at position 2,1 in pUAvec

for soln_space = soln_space22     % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E2, x, v] = RPCalcError2(modifiedVector, u, p, t, N); % Calculate the error at each value
                E22_2 = E2.^2;               % Square the error at each value
                sum_E22_2 = sum(E22_2);      % Sum the squared error at each value
                SSR_22 = [SSR_22 sum_E22_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end

subplot(2,1,2)
title('Likelihood Profile for RP1, parameter UA_amb2');
LLRatio_UAamb2 = 2*log(SSR_22/min(SSR_22));
plot(soln_space22, LLRatio_UAamb2);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
ax = gca;
ax.FontSize = font_size; 
hold off
xlabel('UA_a_m_b_2 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([soln_space22(1) soln_space22(end)]);
x1 = interp1(LLRatio_UAamb2, soln_space22, 0);
zero_point = find(soln_space22 == x1);
x2 = interp1(LLRatio_UAamb2(1:zero_point), soln_space22(1:zero_point), 2.71);
x3 = interp1(LLRatio_UAamb2(zero_point:end), soln_space22(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

% FRIDGE PLANT 3
% For pUA_RP3
%N=373:873;
varyRowIndex    = 1; % Define the row position of the parameter you want to change
varyColumnIndex = 3; % Define the column position of the parameter you want to change

SSR_13 = [];                    % Pre-define a sum of squared residuals
soln_space13 = -6000:1000:30000; % Define the solution space for parameter at position 1,2 in pUAvec

for soln_space = soln_space13     % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E3, x, v] = RPCalcError3(modifiedVector, u, p, t, N); % Calculate the error at each value
                E13_2 = E3.^2;               % Square the error at each value
                sum_E13_2 = sum(E13_2);      % Sum the squared error at each value
                SSR_13 = [SSR_13 sum_E13_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end
font_size = 30;
figure(6)
subplot(2,1,1)
title('Likelihood Profile for RP3, parameter UA_RP3');
LLRatio_UARP3 = 2*log(SSR_13/min(SSR_13));
plot(soln_space13, LLRatio_UARP3);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
ax = gca;
ax.FontSize = font_size; 
hold off
xlabel('UA_R_P_3 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([-500 soln_space13(end)]);
x1 = interp1(LLRatio_UARP3, soln_space13, 0);
zero_point = find(soln_space13 == x1);
x2 = interp1(LLRatio_UARP3(1:zero_point), soln_space13(1:zero_point), 2.71);
x3 = interp1(LLRatio_UARP3(zero_point:end), soln_space13(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

% For pUA_amb3
varyRowIndex    = 2; % Define the row position of the parameter you want to change
varyColumnIndex = 3; % Define the column position of the parameter you want to change

SSR_23 = [];                 % Pre-define a sum of squared residuals
soln_space23 = -4000:2000:18000; % Define the solution space for parameter at position 2,1 in pUAvec

for soln_space = soln_space23     % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E3, x, v] = RPCalcError3(modifiedVector, u, p, t, N); % Calculate the error at each value
                E23_2 = E3.^2;               % Square the error at each value
                sum_E23_2 = sum(E23_2);      % Sum the squared error at each value
                SSR_23 = [SSR_23 sum_E23_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end

subplot(2,1,2)
title('Likelihood Profile for RP3, parameter UA_amb3');
LLRatio_UAamb3 = 2*log(SSR_23/min(SSR_23));
plot(soln_space23, LLRatio_UAamb3);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
hold off
xlabel('UA_a_m_b_3 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([0 500]);
x1 = interp1(LLRatio_UAamb3, soln_space23, 0);
zero_point = find(soln_space23 == x1);
x2 = interp1(LLRatio_UAamb3(1:zero_point), soln_space23(1:zero_point), 2.71);
x3 = interp1(LLRatio_UAamb3(zero_point:end), soln_space23(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

% FRIDGE PLANT 4
% For pUA_RP4

varyRowIndex    = 1; % Define the row position of the parameter you want to change
varyColumnIndex = 4; % Define the column position of the parameter you want to change

SSR_14 = [];                    % Pre-define a sum of squared residuals
soln_space14 = 28000:1000:50000; % Define the solution space for parameter at position 1,2 in pUAvec

for soln_space = soln_space14         % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E4, x, v] = RPCalcError4(modifiedVector, u, p, t, N); % Calculate the error at each value
                E14_2 = E4.^2;               % Square the error at each value
                sum_E14_2 = sum(E14_2);      % Sum the squared error at each value
                SSR_14 = [SSR_14 sum_E14_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end

figure(7)
subplot(2,1,1)
title('Likelihood Profile for RP4, parameter UA_RP4');
LLRatio_UARP4 = 2*log(SSR_14/min(SSR_14));
plot(soln_space14, LLRatio_UARP4);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
ax = gca;
ax.FontSize = font_size; 
hold off
xlabel('UA_R_P_4 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([soln_space14(1) 46000]);
x1 = interp1(LLRatio_UARP4, soln_space14, 0);
zero_point = find(soln_space14 == x1);
x2 = interp1(LLRatio_UARP4(1:zero_point), soln_space14(1:zero_point), 2.71);
x3 = interp1(LLRatio_UARP4(zero_point:end), soln_space14(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

% For pUA_amb4
varyRowIndex    = 2; % Define the row position of the parameter you want to change
varyColumnIndex = 4; % Define the column position of the parameter you want to change

SSR_24 = [];                 % Pre-define a sum of squared residuals
soln_space24 = 0:100:400; % Define the solution space for parameter at position 2,1 in pUAvec

for soln_space = soln_space24         % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E4, x, v] = RPCalcError4(modifiedVector, u, p, t, N); % Calculate the error at each value
                E24_2 = E4.^2;               % Square the error at each value
                sum_E24_2 = sum(E24_2);      % Sum the squared error at each value
                SSR_24 = [SSR_24 sum_E24_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end

subplot(2,1,2)
title('Likelihood Profile for RP4, parameter UA_amb4');
LLRatio_UAamb4 = 2*log(SSR_24/min(SSR_24));
plot(soln_space24, LLRatio_UAamb4);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
hold off
xlabel('UA_a_m_b_4 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([soln_space24(1) soln_space24(end)]);
x1 = interp1(LLRatio_UAamb4, soln_space24, 0);
zero_point = find(soln_space24 == x1);
x2 = interp1(LLRatio_UAamb4(1:zero_point), soln_space24(1:zero_point), 2.71);
x3 = interp1(LLRatio_UAamb4(zero_point:end), soln_space24(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

% FRIDGE PLANT 5
% For pUA_RP5
varyRowIndex    = 1; % Define the row position of the parameter you want to change
varyColumnIndex = 5; % Define the column position of the parameter you want to change

SSR_15 = [];                    % Pre-define a sum of squared residuals
soln_space15 = 24000:100000:900000; % Define the solution space for parameter at position 1,2 in pUAvec

for soln_space = soln_space15         % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E5, x, v] = RPCalcError5(modifiedVector, u, p, t, N); % Calculate the error at each value
                E15_2 = E5.^2;               % Square the error at each value
                sum_E15_2 = sum(E15_2);      % Sum the squared error at each value
                SSR_15 = [SSR_15 sum_E15_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end

figure(8)
subplot(2,1,1)
title('Likelihood Profile for RP5, parameter UA_RP5');
LLRatio_UARP5 = 2*log(SSR_15/min(SSR_15));
plot(soln_space15, LLRatio_UARP5);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
ax = gca;
ax.FontSize = font_size; 
hold off
xlabel('UA_R_P_5 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([soln_space15(1) soln_space15(end)]);
x1 = interp1(LLRatio_UARP5, soln_space15, 0);
zero_point = find(soln_space15 == x1);
x2 = interp1(LLRatio_UARP5(1:zero_point), soln_space15(1:zero_point), 2.71);
x3 = interp1(LLRatio_UARP5(zero_point:end), soln_space15(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

% For pUA_amb5
varyRowIndex    = 2; % Define the row position of the parameter you want to change
varyColumnIndex = 5; % Define the column position of the parameter you want to change

SSR_25 = [];                 % Pre-define a sum of squared residuals
soln_space25 = -2000:500:1000; % Define the solution space for parameter at position 2,1 in pUAvec

for soln_space = soln_space25         % Create a for-loop that loops through different parameter values
    for row = 1:size(pUAvec,1)        % Create a for-loop that loops through different row positions in the vector 
        for column = 1:size(pUAvec,2) % Create a for-loop that loops through different column positions in the vector
            if row == varyRowIndex && column == varyColumnIndex % If the specified vector position is reached
                modifiedVector = pUAvec;                 % Let a new modified vector equal the original vector
                modifiedVector(row,column) = soln_space; % Specified position in the modified vector loops through values in solution space
                [E5, x, v] = RPCalcError5(modifiedVector, u, p, t, N); % Calculate the error at each value
                E25_2 = E5.^2;               % Square the error at each value
                sum_E25_2 = sum(E25_2);      % Sum the squared error at each value
                SSR_25 = [SSR_25 sum_E25_2]; % Build a vector of squared errors for that varied parameter
            else
                pUAvec = pUAvec;          % Otherwise, just let the vector of parameters be equal to the original values
            end
        end
    end
end

subplot(2,1,2)
title('Likelihood Profile for RP5, parameter UA_amb5');
LLRatio_UAamb5 = 2*log(SSR_25/min(SSR_25));
plot(soln_space25, LLRatio_UAamb5);
hold on
yline(2.71,'-',{'Chi-Square Threshold'},'FontSize', font_size);
hold off
xlabel('UA_a_m_b_5 (kJ/Ks)');
ylabel('Negative Log Likelihood Ratio');
xlim([0 500]);
x1 = interp1(LLRatio_UAamb5, soln_space25, 0);
zero_point = find(soln_space25 == x1);
x2 = interp1(LLRatio_UAamb5(1:zero_point), soln_space25(1:zero_point), 2.71);
x3 = interp1(LLRatio_UAamb5(zero_point:end), soln_space25(zero_point:end), 2.71);
xline(x1, '--', {'Optimal Parameter Value'},'FontSize', font_size,'LineWidth', 2);
xline(x2, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
xline(x3, '--', {'90% Confidence Interval'},'FontSize', font_size,'LineWidth', 2);
ax = gca;
ax.FontSize = font_size; 

%% MAPE
% RP1
forecast1 = v.T_RP(:,1);
observed1 = u.T_RP.Values(1:N,1);
ABS1 = abs((forecast1 - observed1)./observed1)*100;
MAPE1 = 1/N*sum(ABS1);

% RP2
forecast2 = v.T_RP(:,2);
observed2 = u.T_RP.Values(1:N,2);
ABS2 = abs((forecast2 - observed2)./observed2)*100;
MAPE2 = 1/N*sum(ABS2);

% RP3
forecast3 = v.T_RP(:,3);
observed3 = u.T_RP.Values(1:N,3);
ABS3 = abs((forecast3 - observed3)./observed3)*100;
MAPE3 = 1/length(N)*sum(ABS3);

% RP4
forecast4 = v.T_RP(:,4);
observed4 = u.T_RP.Values(1:N,4);
ABS4 = abs((forecast4 - observed4)./observed4)*100;
MAPE4 = 1/N*sum(ABS4);

% RP5
forecast5 = v.T_RP(:,5);
observed5 = u.T_RP.Values(1:N,5);
ABS5 = abs((forecast5 - observed5)./observed5)*100;
MAPE5 = 1/N*sum(ABS5);

MAPE = [MAPE1 MAPE2 MAPE3 MAPE4 MAPE5]