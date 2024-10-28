function v = RPIntermediatesKF(x, u, p, t, output, s, response)
t = t(end);
v.m_inRPtot = u.F_inRPtot(t); % kg/s, Total mass flowrate into the fridge plants


% Overall
v.n         = sum(u.s(t),2); % -, Number of RPs in operation calculated by summing each row of s
                              % which is the ON/OFF status of each plant
                              % (either 0 or 1), but is sometimes recorded
                              % as a fraction.

% Calculate mass & volumetric flowrates:
if v.n > 0
    v.F_outRP = (u.F_outCD(t) + u.F_outPT(t))./v.n.*u.s(t); % Note: When the systems are combined, 
                                                            % these two should come from upstream data, 
                                                            % which are already filtered disturbances.
else
    v.F_outRP = u.s(t);
end


v.F_Rec  = [output.MV1(t) output.MV2(t) output.MV3(t) output.MV4(t) output.MV5(t)];


v.F_inRP = v.F_outRP + v.F_Rec;       % L/s,  Volumetric flowrate of the recycle stream for each plant

v.T_inRPtot = ((u.F_outCD(t).*u.T_outCD(t)) + ...
              (u.F_outPT(t).*u.T_outPT(t)))./...
              (u.F_outCD(t)+u.F_outPT(t)); % oC, Temperature of the combined PT and CD outlet streams 
                                           % that join to make the stream entering the fridge plants

v.T_inRP = ((v.F_Rec.*response) + ...
           (v.F_outRP.*v.T_inRPtot))...
           ./ (v.F_Rec + v.F_outRP); % oC, Temperature of the stream entering 
                                     % the evaporator of each individual fridge plant
% WHEN COMBINING THE SYSTEM, REPLACE ALL F_OUTPT AND F_OUTCD WITH THE
% FILTERED DISTURBANCE VARIABLE VALUES FROM UPSTREAM



