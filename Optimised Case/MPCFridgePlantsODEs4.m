function dxdt = MPCFridgePlantsODEs4(p, x, u, t, output, s, i, v)
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
ddt.T4 =     (((v.F_outRP(:,i) + v.F_Rec(:,i)).*x.Tin4)...
            - ((v.F_outRP(:,i) + v.F_Rec(:,i)).*x.T4)...
            + (((p.UA_amb(i) - p.UA_RP(i))./p.C_p).*x.T4))./p.m_RPj...
            + x.Qrefr4(end) - x.Qamb4(end);

ddt.Qrefr4 = 0;
ddt.Qamb4  = 0;
ddt.Tin4   = 0;

dxdt   = S2V(ddt, s.MPCstatefields_RP(:,i));