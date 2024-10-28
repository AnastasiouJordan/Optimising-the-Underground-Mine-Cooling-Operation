function v = RPIntermediates(h_RP, u, p, t)
% Calculate intermediate process variables 
% (i.e. all variables which are neither exogeneous inputs 
%       nor state variables)
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs (measured variables)
%   p: structure of parameters


% Choose whether to use raw data from the plant or simulated data:

% Comment out the first line to use simulated/generated data
% Comment out the second line to use orginal measured plant data

v.m_inRPtot = u.F_inRPtot(t)*p.rho_Water/1000; % kg/s, Total mass flowrate into the fridge plants
%v.m_inRPtot = u.F_inRPtot_generated(t)*p.rho_Water/1000; % kg/s, Total mass flowrate into the fridge plants

% Overall
v.n         = sum(u.s(t),2); % -, Number of RPs in operation calculated by summing each row of s
                              % which is the ON/OFF status of each plant
                              % (either 0 or 1), but is sometimes recorded
                              % as a fraction.

% Calculate mass & volumetric flowrates:

if v.n > 0
    v.m_RP = v.m_inRPtot./v.n.*u.s(t);
else
    v.m_RP = u.s(t);
end

v.m_inRP  = u.F_inRP(t) .* p.rho_Water/1000;  % kg/s, Mass flowrate entering each plant
v.m_outRP = u.F_outRP(t) .* p.rho_Water/1000; % kg/s, Mass flowrate exiting each plant
v.F_RecRP = u.F_outRP(t) - u.F_inRP(t);       % L/s,  Volumetric flowrate of the recycle stream for each plant
v.m_RecRP = v.F_RecRP .* p.rho_Water/1000;    % kg/s, Mass flowrate of the recycle stream for each plant

v.h_inRP = (p.C_p .* (u.T_inRP(t) - p.T_0)) + p.h_0;    % kJ/kg, Enthalpy of the water entering each plant
v.T_RP   = (h_RP' - p.h_0 + (p.C_p .* p.T_0))./p.C_p;   % oC,    Temperature of the water leaving each plant, also measured
v.Q_evap = p.UA_RP .* (v.T_RP - u.T_refr(t)); % kJ/s,  Rate of heat transfer between water in each plant and refrigerant
v.Q_amb  = u.s(t) .* p.UA_amb .* (u.T_amb_meas(t) - v.T_RP); % kJ/s,  Rate of heat transfer between water in each plant and ambient air



