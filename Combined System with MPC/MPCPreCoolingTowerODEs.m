function dxdt = MPCPreCoolingTowerODEs(s, p, x, u, t, output)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.MPCstatefields_PT);
v = PTIntermediatesKF(x, u, p, t, output);

% Calculate state derivatives as structure
ddt.L_PT = ((x.F_AntiSurge(end) + output.MV_SD(t)) - output.MV_PT(t) - v.m_evapPT)./p.m_PTmax*100;
ddt.F_AntiSurge = 0;
ddt.T_AntiSurge = 0;
ddt.T_outSD = 0;
ddt.T_PT = ((output.MV_SD(t).*x.T_outSD(end)) + ...
           (x.F_AntiSurge(end).*x.T_AntiSurge(end)) - ...
           (x.T_PT.*output.MV_SD(t)) - ...
           (x.T_PT.*x.F_AntiSurge(end)) + ...
           ((p.m_Air .* (v.h_inAir - v.h_Air))./p.C_p))./(x.L_PT*p.m_PTmax/100);
ddt.T_amb = 0;
ddt.H_Env = 0;

% Map state derivative structure to vector
dxdt = S2V(ddt,s.MPCstatefields_PT);
