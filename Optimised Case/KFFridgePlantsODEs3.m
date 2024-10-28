 function dxdt = KFFridgePlantsODEs83(s,p,x,u,t,output,i,v)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.KFstatefields_RP(:,i));


% Calculate state derivatives as structure
ddt.T3 =     (((v.F_outRP(:,i) + v.F_Rec(:,i)).*x.Tin3)...
            - ((v.F_outRP(:,i) + v.F_Rec(:,i)).*x.T3)...
            + (((p.UA_amb(i) - p.UA_RP(i))./p.C_p).*x.T3))./p.m_RPj...
            + x.Qrefr3(end) - x.Qamb3(end);

ddt.Qrefr3 = 0;
ddt.Qamb3  = 0;
ddt.Tin3   = 0;


A = [-(v.F_Rec(:,i)./p.m_RPj)-(v.F_outRP(:,i)./p.m_RPj)+(p.UA_amb(i)./(p.C_p*p.m_RPj))-(p.UA_RP(i)./(p.C_p*p.m_RPj)) 1 -1 ((v.F_outRP(:,i)+v.F_Rec(:,i))./p.m_RPj);    
       0 0 0 0;
       0 0 0 0;
       0 0 0 0]; 


ddt.P3 = A.*x.P3.*A' + p.Q3;
ddt.P3 = ddt.P3(1,1);

dxdt   = S2V(ddt, s.KFstatefields_RP(:,i));

