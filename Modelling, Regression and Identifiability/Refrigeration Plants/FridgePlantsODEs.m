 function dxdt = FridgePlantsODEs(p, h_RP, u, t)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables

v = RPIntermediates(h_RP, u, p, t);

% Calculate state derivatives as structure

dxdt =    ((v.m_RP .* v.h_inRP) - (v.m_RP .* h_RP')...
           - v.Q_evap...
           + v.Q_amb)' ./ p.m_RPj;


