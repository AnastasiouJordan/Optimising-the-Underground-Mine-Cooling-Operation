function dxdt = KFPreCoolingTowerODEs(s, p, x, u, t, output)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.KFstatefields_PT);
v = PTIntermediatesKF(x, u, p, t, output);

% Calculate state derivatives as structure
ddt.L_PT = ((x.F_AntiSurge(end) + output.MV_SD(t)) - output.MV_PT(t) - v.m_evapPT)./p.m_PTmax*100;
ddt.F_AntiSurge = 0;
ddt.T_AntiSurge = 0;
ddt.T_outSD = 0;
ddt.T_PT = ((output.MV_SD(t).*x.T_outSD(end)) + ...
           (x.F_AntiSurge(end).*x.T_AntiSurge(end)) - ...
           (x.T_PT(end).*output.MV_SD(t)) - ...
           (x.T_PT(end).*x.F_AntiSurge(end)) + ...
           ((p.m_Air .* (v.h_inAir - v.h_Air))./p.C_p))./(x.L_PT*p.m_PTmax/100);
ddt.T_amb = 0;
ddt.H_Env = 0;

A = [0 p.c_L_PT 0 0 0 0 0;    
     (-100/(x.L_PT.^2.*(x.L_PT*p.m_PTmax/100))).*((output.MV_SD(t).*x.T_outSD(end))+(x.F_AntiSurge(end).*x.T_AntiSurge(end))-(x.T_PT.*output.MV_SD(t))-(x.T_PT.*x.F_AntiSurge(end))+((p.m_Air.*(v.h_inAir - v.h_Air))./p.C_p)) ((x.T_AntiSurge(end)-x.T_PT)./(x.L_PT*p.m_PTmax/100)) (x.F_AntiSurge(end)./(x.L_PT*p.m_PTmax/100)) (output.MV_SD(t)./(x.L_PT*p.m_PTmax/100)) (-(output.MV_SD(t)+x.F_AntiSurge(end))./(x.L_PT*p.m_PTmax/100)) 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 0 0 0 0];     


ddt.P_PT = A.*x.P_PT.*A' + p.Q_PT;
ddt.P_PT = ddt.P_PT(1,1);

dxdt = S2V(ddt,s.KFstatefields_PT);

