function dxdt = PreCoolingTowersODEs(s, p, x, u, t, output, z)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   x_vec: vector of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters
%   s: structure of state field names

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.statefields);
v = PTIntermediates(x, u, p, t);

% Calculate state derivatives as structure
ddt.L = (u.F_in_generated(t) - output.MV(t) - v.m_evapPT)./p.m_PTmax*100;
ddt.h_PT = ((u.F_in_generated(t).*(v.h_inPT - x.h_PT')) + ...
           (p.m_Air .* (v.h_inAir - v.h_Air)))./ (x.L/100*p.m_PTmax); % Energy balance
                                                         % for the water
                                                         % in the PTs


% Map state derivative structure to vector
dxdt = S2V(ddt, s.statefields);
