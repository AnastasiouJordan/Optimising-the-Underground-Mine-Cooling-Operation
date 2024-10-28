function dxdt = MPCChillDamODEs(s, p, x, u, t, output,v)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.MPCstatefields_CD);

% Calculate mass & volumetric flowrates:
v.m_1       = v.F_outRPtot(end) - x.F_UG(end); % kg/s, Inlet mass flowrate based on filtered measurements
v.m_2       = output.MV_CD(t);            % kg/s, Outlet mass flowrate based on filtered measurements
v.m_3       = x.F_Ice(end);               % kg/s, Outlet Ice plant flowrate based on filtered steps         

v.T_1 = v.T_outRPtot(end).*(v.m_1 > 0) + x.T_CD'.*(v.m_1 < 0); % If the flowrate into the CD is positive, let the temp be equal to that of the water leaving the Fridge Plants.
                                                   % If the flowrate into the CD is negative, let the temp be equal to that of the water in the Chill Dam.
v.T_2 = x.T_CD'.*(v.m_2 > 0) + v.T_PT.*(v.m_2 < 0); % If the flowrate out of the CD is positive, let the temp be equal to that of the water in the Chill Dam.
                                                   % If the flowrate out of the CD is negative, let the temp be equal to that of the water leaving the Pre-Cooling Towers.
v.T_3 = x.T_CD';

% Calculate state derivatives as structure
ddt.L_CD = (v.m_1 - v.m_2 - v.m_3)./p.m_CDmax*100;

% ddt.T_CD = (v.T_1.*v.m_1 - (v.m_1 .* x.T_CD') - ...
%             v.T_2.*v.m_2 + (v.m_2 .* x.T_CD') - ...
%             v.T_3.*v.m_3 + (v.m_3 .* x.T_CD') + ...
%             (p.UA_CD .* (v.T_amb - x.T_CD'))./p.C_p) ./ (x.L_CD'./100.*p.m_CDmax);

ddt.T_CD = (v.T_1.*v.m_1 - ...
            v.T_2.*v.m_2 - ...
            v.T_3.*v.m_3 + ...
            (p.UA_CD .* (v.T_amb - x.T_CD))./p.C_p) ./ (x.L_CD./100.*p.m_CDmax);

ddt.F_Ice = 0;
ddt.F_UG  = 0;

% Map state derivative structure to vector
dxdt = S2V(ddt,s.MPCstatefields_CD);