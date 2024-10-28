function dxdt = MPCPreCoolingTowerODEs(s, p, x, u, t, output)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.MPCstatefields);
v = PTIntermediates(x, u, p, t);

% Calculate state derivatives as structure
ddt.L = (x.F_in(end) - output.MV(t) - v.m_evapPT)./p.m_PTmax*100;
ddt.F_in = 0;
% ddt.h_PT = ((x.F_in(end).*(v.h_inPT - x.h_PT')) + ...
%            (p.m_Air .* (v.h_inAir - v.h_Air)))./ (x.L/100*p.m_PTmax); % Energy balance
%                                                          % for the water
%                                                          % in the PTs
ddt.h_PT = 0;

% Map state derivative structure to vector
dxdt = S2V(ddt,s.MPCstatefields);
