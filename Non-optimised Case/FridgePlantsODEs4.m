function dxdt = FridgePlantsODEs4(s,p,x,u,t,output,i,response)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

x = V2S(x, s.statefields_RP(:,i));

% Calculate the intermediates
v.m_inRPtot = output.MV_PT(t) + output.MV_CD(t); % kg/s, Total mass flowrate into the fridge plants
v.n         = output.MV_Econ(end);  % -, Number of RPs in operation calculated by summing each row of s
                              % which is the ON/OFF status of each plant
                              % (either 0 or 1), but is sometimes recorded
                              % as a fraction.
u_s = [1*(v.n>=1) 1*(v.n>1) 1*(v.n>2) 1*(v.n>3) 1*(v.n>4)];

v.n = sum(u_s,2);

if v.n > 0
    v.F_outRP = (output.MV_CD(t) + output.MV_PT(t))./v.n.*u_s; % Note: When the systems are combined, 
                                                            % these two should come from upstream data, 
                                                            % which are already filtered disturbances.
else
    v.F_outRP = u_s+0.001;
end

v.F_Rec     = output.MV1(t);
v.F_inRP    = v.F_outRP(:,i) + v.F_Rec;       % L/s,  Volumetric flowrate of the recycle stream for each plant
v.T_inRPtot = ((output.MV_CD(t).*response.CD(2,end)') + ...
              (output.MV_PT(t).*response.PT(2,end)'))./...
              (output.MV_CD(t) + output.MV_PT(t)); % oC, Temperature of the combined PT and CD outlet streams 
                                           % that join to make the stream entering the fridge plants
v.T_inRP    = ((v.F_Rec.*response.RP(end,i)) + ...
                (v.F_outRP(:,i).*v.T_inRPtot))...
                ./ (v.F_Rec + v.F_outRP(:,i)); % oC, Temperature of the stream entering 
                                                  % the evaporator of each individual fridge plant


% WHEN COMBINING THE SYSTEM, REPLACE ALL F_OUTPT AND F_OUTCD WITH THE
% FILTERED DISTURBANCE VARIABLE VALUES FROM UPSTREAM

v.Q_refr_generated = u.Q_refr_generated(t);
v.Q_amb_generated  = u.Q_amb_generated(t);
v.T_in_generated   = u.T_in_generatedRP(t);

% Calculate state derivatives as structure
ddt.T4 =     (((v.F_outRP(:,i) + v.F_Rec).*v.T_in_generated(:,i))...
            - ((v.F_outRP(:,i) + v.F_Rec).*x.T4)...
            + (((p.UA_amb(i) - p.UA_RP(i))./p.C_p).*x.T4))./p.m_RPj...
            + v.Q_refr_generated(:,i) - v.Q_amb_generated(:,i);

dxdt   = S2V(ddt, s.statefields_RP(:,i));

