function dxdt = PreCoolingTowersODEs(s, p, x, u, t, output, z)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   x_vec: vector of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters
%   s: structure of state field names

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.statefields_PT);
v = PTIntermediates(x, u, p, t, output);

% Calculate state derivatives as structure
ddt.L_PT = ((u.F_AntiSurge_generated(t) + output.MV_SD(t)) - output.MV_PT(t) - v.m_evapPT)./p.m_PTmax*100;
ddt.T_PT = ((output.MV_SD(t).*u.T_out_generatedSD(t)) + ...
           (u.F_AntiSurge_generated(t).*u.T_AntiSurge_generated(t)) - ...
           (x.T_PT.*output.MV_SD(t)) - ...
           (x.T_PT.*u.F_AntiSurge_generated(t)) + ...
           ((p.m_Air .* (v.h_inAir - v.h_Air))./p.C_p))./(x.L_PT*p.m_PTmax/100);


% Map state derivative structure to vector
dxdt = S2V(ddt, s.statefields_PT);
