function dxdt = SquareDamODEs(s, p, x_vec, u, t)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x_vec, s.statefields);
v = SDIntermediates(x, u, p, t);

% Calculate state derivatives as structure
ddt.m_SD = u.F_in_filtered(t) - u.F_out_filtered(t) - p.m_evapSD;


% Map state derivative structure to vector
dxdt = S2V(ddt,s.statefields);
