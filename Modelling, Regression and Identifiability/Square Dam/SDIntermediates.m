function v = SDIntermediates(x, u, p, t)
% Calculate intermediate process variables 
% (i.e. all variables which are neither exogeneous inputs 
%       nor state variables)
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs (measured variables)
%   p: structure of parameters


v.m_inSD  = u.F_inSD(t) * p.rho_Water / 1000;  % kg/s, Mass flowrate into the SD
v.m_outSD = u.F_outSD(t) * p.rho_Water / 1000; % kg/s, Mass flowrate out of the SD

v.L_SD     = x.m_SD' .* 100 / p.m_SDmax;         % %, Level in the SD, also measured