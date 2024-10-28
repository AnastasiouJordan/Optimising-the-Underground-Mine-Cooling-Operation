function v = RPIntermediatesKF(x, u, p, t, output, s, response)

% This function is going to be used in the KF, using only state estimates.

t = t(end);
v.m_inRPtot = output.MV_PT(t) + output.MV_CD(t); % kg/s, Total mass flowrate into the fridge plants


% Overall
v.n         = output.MV_Econ(end); % -, Number of RPs in operation calculated by summing each row of s
                              % which is the ON/OFF status of each plant
                              % (either 0 or 1), but is sometimes recorded
                              % as a fraction.
 
u_s = [1*(v.n>=1) 1*(v.n>1) 1*(v.n>2) 1*(v.n>3) 1*(v.n>4)];

v.n = sum(u_s,2);

% Calculate mass & volumetric flowrates:
if v.n > 0
    v.F_outRP = (output.MV_CD(t) + output.MV_PT(t))./v.n.*u_s; % Note: When the systems are combined, 
                                                            % these two should come from upstream data, 
                                                            % which are already filtered disturbances.
else
    v.F_outRP = u_s+0.001;
end


v.F_Rec  = [output.MV1(t) output.MV2(t) output.MV3(t) output.MV4(t) output.MV5(t)];


v.F_inRP = v.F_outRP + v.F_Rec;       % L/s,  Volumetric flowrate of the recycle stream for each plant

v.T_inRPtot = ((output.MV_CD(t).*x.T_CD(end)') + ...
              (output.MV_PT(t).*x.T_PT(end)'))./...
              (output.MV_CD(t) + output.MV_PT(t)); % oC, Temperature of the combined PT and CD outlet streams 
                                           % that join to make the stream entering the fridge plants


v.T_RP   = [x.T1(end) x.T2(end) x.T3(end) x.T4(end) x.T5(end)];

v.T_inRP = ((v.F_Rec.*v.T_RP) + ...
           (v.F_outRP.*v.T_inRPtot))...
           ./ (v.F_Rec + v.F_outRP); % oC, Temperature of the stream entering 
                                     % the evaporator of each individual fridge plant                             
                                     

v.F_outRPtot = (v.F_outRP(:,1) + v.F_outRP(:,2) + v.F_outRP(:,3) + v.F_outRP(:,4) + v.F_outRP(:,5));
v.T_outRPtot = ((v.F_outRP(:,1).*x.T1) + (v.F_outRP(:,2).*x.T2) + (v.F_outRP(:,3).*x.T3) + (v.F_outRP(:,4).*x.T4)...
    + (v.F_outRP(:,5).*x.T5)) ./ v.F_outRPtot;
