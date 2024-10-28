function dxdt = MPCChillDamODEs(s, p, x, u, t, output)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.MPCstatefields);
v = CDIntermediates(x, u, p, t, output);

% Calculate mass & volumetric flowrates:
v.m_1       = x.F_in(end);          % kg/s, Inlet mass flowrate based on filtered measurements
v.m_2       = output.MV(t);       % kg/s, Outlet mass flowrate based on filtered measurements
v.m_3       = u.F_Ice_generated(t); % kg/s, Outlet Ice plant flowrate based on filtered steps         

v.H_1 = v.h_RP.*(v.m_1 > 0) + x.h_CD'.*(v.m_1 < 0); % If the flowrate into the CD is positive, let the enthalpy be equal to that of the water leaving the Fridge Plants.
                                                   % If the flowrate into the CD is negative, let the enthalpy be equal to that of the water in the Chill Dam.
v.H_2 = x.h_CD'.*(v.m_2 > 0) + v.h_PT.*(v.m_2 < 0); % If the flowrate out of the CD is positive, let the enthalpy be equal to that of the water in the Chill Dam.
                                                   % If the flowrate out of the CD is negative, let the enthalpy be equal to that of the water leaving the Pre-Cooling Towers.
v.H_3 = x.h_CD'.*v.m_3;

% Calculate state derivatives as structure
%ddt.m_CD = u.F_out_PT_filtered(t) - u.F_UG_filtered(t) - u.F_Ice_filtered(t);
ddt.L = (v.m_1 - v.m_2 - v.m_3)./p.m_CDmax*100;
ddt.F_in = 0;
% ddt.h_CD = (v.H_1 - (v.m_1 .* x.h_CD) - ...
%             v.H_2 + (v.m_2 .* x.h_CD) - ...
%             v.H_3 + (v.m_3 .* x.h_CD) + ...
%             (p.UA_CD .* (u.T_Env(t) - v.T_CD))) ./ (x.L/100*p.m_CDmax);

% ddt.h_CD = (v.H_1 - (v.m_1 .* x.h_CD) - ...
%             v.H_2 + (v.m_2 .* x.h_CD) + ...  
%             (p.UA_CD .* (u.T_Env(t) - v.T_CD))) ./(x.L/100*p.m_CDmax);
ddt.h_CD = 0;

% Map state derivative structure to vector
dxdt = S2V(ddt, s.MPCstatefields);