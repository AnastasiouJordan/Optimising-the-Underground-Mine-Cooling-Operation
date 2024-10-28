function dxdt = ChillDamODEs(s, p, x, u, t, output,response)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.statefields_CD);
v = RPIntermediates(x, u, p, t, output, s, response);

% Calculate mass & volumetric flowrates:
v.m_1       = v.F_outRPtot(end) - u.F_UG_generated(t);  % kg/s, Inlet mass flowrate based on filtered measurements
v.m_2       = output.MV_CD(t);            % kg/s, Outlet mass flowrate based on filtered measurements
v.m_3       = u.F_Ice_generated(t); % kg/s, Outlet Ice plant flowrate based on filtered steps         

v.T_1 = v.T_outRPtot(end).*(v.m_1 > 0) + x.T_CD'.*(v.m_1 < 0); % If the flowrate into the CD is positive, let the temp be equal to that of the water leaving the Fridge Plants.
                                                   % If the flowrate into the CD is negative, let the temp be equal to that of the water in the Chill Dam.
v.T_2 = x.T_CD'.*(v.m_2 > 0) + response.PT(2,end).*(v.m_2 < 0); % If the flowrate out of the CD is positive, let the temp be equal to that of the water in the Chill Dam.
                                                   % If the flowrate out of the CD is negative, let the temp be equal to that of the water leaving the Pre-Cooling Towers.
v.T_3 = x.T_CD';

% Calculate state derivatives as structure
ddt.L_CD = (v.m_1 - v.m_2 - v.m_3)./p.m_CDmax*100;

% ddt.T_CD = (v.T_1.*v.m_1 - (v.m_1 .* x.T_CD') - ...
%             v.T_2.*v.m_2 + (v.m_2 .* x.T_CD') - ...
%             v.T_3.*v.m_3 + (v.m_3 .* x.T_CD') + ...
%             (p.UA_CD .* (u.T_amb_generated(t) - x.T_CD'))./p.C_p) ./ (x.L_CD'/100*p.m_CDmax);

ddt.T_CD = (v.T_1.*v.m_1 - ...
            v.T_2.*v.m_2 - ...
            v.T_3.*v.m_3 + ...
            (p.UA_CD .* (u.T_amb_generated(t) - x.T_CD'))./p.C_p) ./ (x.L_CD'/100*p.m_CDmax);

% Map state derivative structure to vector
dxdt = S2V(ddt, s.statefields_CD);

