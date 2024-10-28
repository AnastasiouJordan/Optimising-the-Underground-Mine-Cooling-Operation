function dxdt = MPCFridgePlantsODEs1(p, x, u, t, output, s, i, v)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.MPCstatefields_RP(:,i));


% Calculate state derivatives as structure
ddt.T1 =     (((v.F_outRP(:,i) + v.F_Rec(:,i)).*x.Tin1)...
            - ((v.F_outRP(:,i) + v.F_Rec(:,i)).*x.T1)...
            + (((p.UA_amb(i) - p.UA_RP(i))./p.C_p).*x.T1))./p.m_RPj...
            + x.Qrefr1(end) - x.Qamb1(end);

ddt.Qrefr1 = 0;
ddt.Qamb1  = 0;
ddt.Tin1   = 0;


dxdt   = S2V(ddt, s.MPCstatefields_RP(:,i));


