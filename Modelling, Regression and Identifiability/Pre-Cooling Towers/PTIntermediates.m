function v = PTIntermediates(x, u, p, t)
% Calculate intermediate process variables 
% (i.e. all variables which are neither exogeneous inputs 
%       nor state variables)
%
% The function requires the following process variables as inputs:
%   x: structure of state variables
%   u: structure of exogeneous inputs (measured variables)
%   p: structure of parameters
%   t: time


% Calculate mass & volumetric flowrates:
v.F_inPT  = u.F_in_filtered(t);                     % L/s,  Inlet volumetric flowrate filtered measurements
v.m_inPT  = v.F_inPT * p.rho_Water/1000;            % kg/s, Inlet mass flowrate based on filtered measurements
v.F_outPT = u.F_out_filtered(t);                    % L/s,  Outlet volumetric flowrate filtered measurements
v.m_outPT = v.F_outPT * p.rho_Water/1000; % kg/s, Outlet mass flowrate based on filtered measurements

% Calculate enthalpy for the air entering at Env conditions
v.T_Env_Kel = u.T_Env(t) + 273.15; % K, Convert measured dry bulb 
                                   % environmental temperature
                                   % from oC to Kelvin for use in the
                                   % P_wsTEnv equation

v.P_wsTEnv  = exp((p.C1./v.T_Env_Kel) + p.C2 + (p.C3*v.T_Env_Kel) + ...
              (p.C4*v.T_Env_Kel.^2) + (p.C5*v.T_Env_Kel.^3) + ...
              (p.C6*log(v.T_Env_Kel)))./1000; % kPa, Vapour saturation pressure
                                              % of condensed water entering PTs                                              


v.P_wEnv    = v.P_wsTEnv.*u.H_Env(t)./100; % kPa, Partial pressure of water
                                           % vapour entering PTs


v.W_Env     = (p.M_r*v.P_wEnv)./(p.P_Tot - v.P_wEnv); % kg/kg, Mixing ratio

v.h_inAir   = u.T_Env(t).*(p.A + (p.B.*v.W_Env))...
              + (p.C.*v.W_Env); % kJ/kg, Enthalpy of air entering PTs


% Calculate enthalpy for air leaving at PT conditions
v.T_PT     = (x.h_PT' - p.h_0 + (p.C_p * p.T_0)) ./ p.C_p; % oC, Temperature of 
                                                           % water leaving PTs.
                                                           % Also measured as
                                                           % u.T_PT.


v.T_PT_Kel = v.T_PT + 273.15; % Convert calculated outlet PT temperature
                              % of water from oC to Kelvin for use in the
                              % P_wsTPT equation.

% v.P_wsTPT = exp((p.C1./v.T_PT_Kel) + p.C2 + (p.C3*v.T_PT_Kel) + ...
%                    (p.C4*v.T_PT_Kel.^2) + (p.C5*v.T_PT_Kel.^3) + ...
%                    (p.C6*log(v.T_PT_Kel)))/1000; % kPa, Vapour saturation pressure
%                                                  % of condensed water leaving
v.T_PT_Kel = max(273, v.T_PT_Kel);
% v.P_wsTPT = exp((p.C1./(max(5, v.T_PT_Kel))) + p.C2 + (p.C3*v.T_PT_Kel) + ...
%              (p.C4*v.T_PT_Kel.^2) + (p.C5*v.T_PT_Kel.^3) + ...
%              (p.C6*log(max(5, v.T_PT_Kel))))/1000;
v.P_wsTPT = exp((p.C1./v.T_PT_Kel) + p.C2 + (p.C3*v.T_PT_Kel) + ...
             (p.C4*v.T_PT_Kel.^2) + (p.C5*v.T_PT_Kel.^3) + ...
             (p.C6*log(v.T_PT_Kel)))/1000;
%v.P_wsTPT = 1.5;


v.P_wPT    = v.P_wsTPT.*p.H_outAir./100; % kPa, Partial pressure of water
                                        % vapour leaving PTs.

v.W_PT     = (p.M_r*v.P_wPT)./(p.P_Tot - v.P_wPT); % kg/kg, Mixing ratio.

v.h_Air    = v.T_PT .*(p.A + (p.B.*v.W_PT))...
             + (p.C.*v.W_PT); % kJ/kg, Enthalpy of air leaving PTs.
                              % Here the temperature in oC is used.

% Evaporation
v.m_evapPT = (p.H_outAir/100).*(p.m_Air*0.00216679*v.P_wsTPT)...
             ./(p.rho_Air*v.T_PT_Kel); % kg/s, Mass flowrate of water leaving
                                       % the PTs as evaporation. This
                                       % equation is derived by combining
                                       % absolute and relative humidity
                                       % equations, where the constant
                                       % '0.00216679' comes from the
                                       % relationship between absolute and
                                       % relative humidity and originally
                                       % is in units of kgK/J. It is then
                                       % multiplied by 1000 during unit
                                       % conversions of pressure and
                                       % flowrates, and is not used in any
                                       % other equations. Dividing the
                                       % humidity by 100 is simply a
                                       % conversion from percentage to
                                       % decimal.

% Water entering PTs
v.T_inPT  = (u.T_outSD(t).*u.F_outSD(t) + u.T_AntiSurge(t)...
            .*u.F_AntiSurge(t))./v.F_inPT; % oC, Temperature of the water entering
                                 % the PTs. Calculated based on the
                                 % measured temperatures of the square dam
                                 % outlet and the antisurge flow which,
                                 % together, form the inlet to the PTs.
                                 % This was where division by zero was
                                 % occuring. Denominator changed from 
                                 % u.F_outSD(t) + u.F_AntiSurge(t) to
                                 % v.F_inPT.

v.h_inPT  = p.C_p*(v.T_inPT - p.T_0) + p.h_0; % kJ/kg, Enthalpy of the
                                              % water entering the PTs.
                                              % Based on the inlet
                                              % temperature and reference
                                              % enthalpy and temperature
                                              % values for T = 0 - 60 oC.


% Calculate levels:
v.L_PTA = u.L_PT_filtered(t)*100./p.height_PT; % %, Conversion of Kalman Filter level estimates
v.L_PTB = u.L_PT_filtered(t)*100./p.height_PT; % %, Conversion of Kalman Filter level estimates


v.m_PT  = (u.L_PTA(t)./100*p.m_PTAmax) + (u.L_PTA(t)./100*p.m_PTBmax);


v.L_PT  = x.m_PT'./(p.m_PTAmax + p.m_PTBmax).*100;