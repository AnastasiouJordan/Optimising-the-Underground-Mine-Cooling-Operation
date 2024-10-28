function dxdt = PreCoolingTowersODEs(s, p, x_vec, u, t)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   x_vec: vector of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters
%   s: structure of state field names

% Map state vector to structure and calculate intermediate variables

x = V2S(x_vec, s.statefields);
v = PTIntermediates(x, u, p, t);


% Calculate state derivatives as structure

ddt.m_PT = u.F_in_filtered(t) - u.F_out_filtered(t) - v.m_evapPT;
ddt.h_PT = ((v.m_inPT.*(v.h_inPT - x.h_PT')) + ...
           (p.m_Air .* (v.h_inAir - v.h_Air)))./ x.m_PT'; % Energy balance
                                                         % for the water
                                                         % in the PTs


% Map state derivative structure to vector
dxdt = S2V(ddt, s.statefields);
%dxdt = dxdt(:,1);

%dbstop in PreCoolingTowersODEs at 13 if isnan(x.h_PT)