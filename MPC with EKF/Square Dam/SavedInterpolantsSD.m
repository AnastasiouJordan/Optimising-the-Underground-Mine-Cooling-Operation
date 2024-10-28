clc
clear

%% Read in data
SD_meas = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','SD','Range','C11670:L62756');
T = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','CD','Range','B11670:Q62756');
t = T{:,1};
t = t - 700020;

%% Define variables from data
% Required
L_SD    =  SD_meas{:,6}; % %,   Level in the SD (also calculated)
F_inSD  =  SD_meas{:,5}; % L/s, Inlet volmuetric flowrate to the SD
F_outSD =  SD_meas{:,8}; % L/s, Outlet volumetric flowrate from the CD

% Redundant
T_outSD = SD_meas{:,9};  % oC, Temperature of the water leaving the SD
U_outSD = SD_meas{:,10}; % %,  Valve position on SD outlet

%  for i = 1:51087
 
% if L_SD(i,1) > quantile(L_SD,[0.75]);
%     L_SD(i,1) = quantile(L_SD,[0.75]);
% end
% if L_SD(i,1) < quantile(L_SD,[0.25]);
%     L_SD(i,1) = quantile(L_SD,[0.25]);
% end
% 
%  if F_inSD(i,1) > quantile(F_inSD,[0.75]);
%      F_inSD(i,1) = quantile(F_inSD,[0.75]);
%  end
%  if F_inSD(i,1) < quantile(F_inSD,[0.25]);
%      F_inSD(i,1) = quantile(F_inSD,[0.25]);
%  end
% % 
%  if F_outSD(i,1) > quantile(F_outSD,[0.75]);
%      F_outSD(i,1) = quantile(F_outSD,[0.75]);
%  end
%  if F_outSD(i,1) < quantile(F_outSD,[0.25]);
%      F_outSD(i,1) = quantile(F_outSD,[0.25]);
%  end
%  
%  end
%% Stochastically model the evaporation from the SD
N = 51087;
m_evapSD = zeros(N, 1); % Just to define the variable and ensure that
                        % the first value in the for-loop is taken as 0.
a   = 0.99; % Indication of how closely correlated the next value is with
            % respect to the previous value
b   = 0.1;  % Indication of the weighting given to the randomness added
            % in each step
% Note: ambient wind speeds are measured, however they are measured in m/s.
% Because this is such a small quantity, we will model it stochastically.
% Assuming wind speeds of between 100 - 250 kg/s:
% min_evap = 0.00001; % kg/s, Rate of water evaporation at 5oC and 1% humidity
% max_evap = 0.008;   % kg/s, Rate of water evaporation at 35oC and 100% humidity
% Therefore specifiy the range of random values to be generated:

max_r = 0.0001;
min_r = 0.00001;

for i = 2:N
    random(i) = (max_r - min_r).*rand(1,1) + min_r;
    m_evapSD(i) = m_evapSD(i-1)*a + random(i)*b; % The next value of m_evapSD should be
end                                               % based on the previous value with
                                                  % some random error added to it
% Now produced random values for m_evapSD which 'follow on' from one another,
% we must ensure that they are non-negative and more typical of an evaporation flowrate.
% We therefore add a mean value as follows:

m_evapSD = m_evapSD + 0.001; % Here, the mean value of the water evaporation becomes 0.001, and
               % all the simulated values in the for-loop contribute
               % towards noisy measurements or errors around the mean.

% In the next section, we can now define m_evapSD as an exogenous input 
% structure and make it continuous with time, as it is currently a set of
% points.          

%% Create continuous data
u.L_SD    = griddedInterpolant(t, L_SD);
u.F_inSD  = griddedInterpolant(t, F_inSD);
u.F_outSD = griddedInterpolant(t, F_outSD);
u.T_outSD = griddedInterpolant(t, T_outSD);
u.U_outSD = griddedInterpolant(t, U_outSD);
u.m_evapSD = griddedInterpolant(t, m_evapSD);

%% Define field names for exogenous data
n.exogenousfields = {'L_SD', 'F_inSD', 'F_outSD',...
                     'T_outSD', 'U_outSD', 'm_evapSD'};

%% Save the data to be used in the main file
save SavedInterpolantsSD.mat u n t

clear all