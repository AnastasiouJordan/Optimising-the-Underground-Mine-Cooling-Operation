clc
clear

CD_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','CD','Range','C11670:L62756');
RP_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','RP','Range','C11672:CG62758');
PT_meas   = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','PT','Range','C11670:Q62756');
ENV_meas  = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','Env','Range','C11670:Q62756');

% Overall
F_inRP   =  RP_meas{:,5};  % L/s, Total inlet volumetric flowrate to RPs
F_outPT  =  PT_meas{:,10}; % L/s, Outlet volmuetric flowrate from PTs
F_outCD  =  CD_meas{:,8};  % L/s, Outlet volumetric flowrate from CD
T_amb_meas = ENV_meas{:,4}; % oC, Ambient temperature measured on site


T = readtable('2022.05.06 - 2022.06.14.xlsx','Sheet','CD','Range','B11670:Q62756');
t = T{:,1};
t = t - 700020;

% Model the ambient temperature surrounding the RPs stochastically
N = 51087;
T_amb = zeros(N, 1); % Just to define the variable and ensure that
                     % the first value in the for-loop is taken as 0.
a   = 0.9; % Indication of how closely correlated the next value is with
            % respect to the previous value
b   = 0.1;  % Indication of the weighting given to the randomness added
            % in each step

for i = 2:N
    T_amb(i) = T_amb(i-1)*a + randn*b; % The next value of T_amb should be
end                                               % based on the previous value with
                                                  % some random error added to it
% Now produced random values for T_amb which 'follow on' from one another,
% we must ensure that they are non-negative and more typical of an ambient temp.
% We therefore add a mean value as follows:

T_amb = T_amb + 20; % Here, the mean value of the water evaporation becomes 20, and
                    % all the simulated values in the for-loop contribute
                    % towards noisy measurements or errors around the mean.


u.F_inRPtot = griddedInterpolant(t, F_inRP);
u.F_outPT   = griddedInterpolant(t, F_outPT); 
u.F_outCD   = griddedInterpolant(t, F_outCD);
u.T_amb     = griddedInterpolant(t, T_amb);
u.T_amb_meas = griddedInterpolant(t, T_amb_meas);

% Per fridge plant, j
% On/off statuses
for j = 1:5
    k = (15*j) + 1;
    sj(:,j) = RP_meas{:, k};
end
sj(sj <= 0.1) = 0;
sj(sj > 0.1)  = 1;
u.s = griddedInterpolant(t, sj); 
plot(t, u.s.Values(:,1), t, sj(:,1)) 
legend('u.s', 'sj')


% Inlet volumetric flowrates
for j = 1:5
    k = (15*j) - 2;
    F_inRPj(:,j) = RP_meas{:, k};
end
u.F_inRP = griddedInterpolant(t, F_inRPj); 


% Outlet volumetric flowrates
for j = 1:5
    k = (15*j);
    F_outRPj(:,j) = RP_meas{:, k}; 
end
u.F_outRP = griddedInterpolant(t, F_outRPj);

% Inlet temperatures (oC)
for j = 1:5
    k = (15*j) - 4;
    T_inRPj(:,j) = RP_meas{:, k};
end
u.T_inRP = griddedInterpolant(t, T_inRPj); 

% Outlet temperatures
for j = 1:5
    k = (15*j) - 1;
    T_RPj(:,j) = RP_meas{:, k};
end
u.T_RP = griddedInterpolant(t, T_RPj); 


% Refrigerant temperatures (liquid)
for j = 1:5
    k = (15*j) + 8;
    T_refrj(:,j) = RP_meas{:, k};
end
u.T_refr = griddedInterpolant(t, T_refrj); 


n.exogenousfields = {'F_inRP', 'F_outPT', 'F_outCD', 'T_amb',...
                     'sj', 'F_inRPj', 'F_outRPj','T_inRPj'...
                     'T_RPj', 'T_refrj'};

save SavedInterpolantsRP.mat u n t
%clear all



