function dxdt = MPCFridgePlantsODEs1(p, x, u, t, output, s, i)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.MPCstatefields(:,i));

% Calculate the intermediates
v.m_inRPtot = u.F_inRPtot(t); % kg/s, Total mass flowrate into the fridge plants
v.n         = sum(u.s(t),2);  % -, Number of RPs in operation calculated by summing each row of s
                              % which is the ON/OFF status of each plant
                              % (either 0 or 1), but is sometimes recorded
                              % as a fraction.
if v.n > 0
    v.F_outRP = (u.F_outCD(t) + u.F_outPT(t))./v.n.*u.s(t); % Note: When the systems are combined, 
                                                            % these two should come from upstream data, 
                                                            % which are already filtered disturbances.
else
    v.F_outRP = u.s(t)+0.001;
end

v.F_Rec     = output.MV1(t);
v.F_inRP    = v.F_outRP(:,i) + v.F_Rec;       % L/s,  Volumetric flowrate of the recycle stream for each plant
v.T_inRPtot = ((u.F_outCD(t).*u.T_outCD(t)) + ...
              (u.F_outPT(t).*u.T_outPT(t)))./...
              (u.F_outCD(t)+u.F_outPT(t)); % oC, Temperature of the combined PT and CD outlet streams 
                                           % that join to make the stream entering the fridge plants
v.T_inRP    = ((v.F_Rec.*x.T1) + ...
              (v.F_outRP(:,i).*v.T_inRPtot))...
              ./ (v.F_Rec + v.F_outRP(:,i)); % oC, Temperature of the stream entering 
                                                  % the evaporator of each individual fridge plant
% WHEN COMBINING THE SYSTEM, REPLACE ALL F_OUTPT AND F_OUTCD WITH THE
% FILTERED DISTURBANCE VARIABLE VALUES FROM UPSTREAM


% Calculate state derivatives as structure
ddt.T1 =     (((v.F_outRP(:,i) + v.F_Rec).*x.Tin1)...
            - ((v.F_outRP(:,i) + v.F_Rec).*x.T1)...
            + (((p.UA_amb(i) - p.UA_RP(i))./p.C_p).*x.T1))./p.m_RPj...
            + x.Qrefr1(end) - x.Qamb1(end);

ddt.Qrefr1 = 0;
ddt.Qamb1  = 0;
ddt.Tin1   = 0;


dxdt   = S2V(ddt, s.MPCstatefields(:,i));


