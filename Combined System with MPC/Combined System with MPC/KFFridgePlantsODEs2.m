 function dxdt = KFFridgePlantsODEs2(s,p,x,u,t,output,i,v)
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
ddt.T2 =     (((v.F_outRP(:,i) + v.F_Rec(:,i)).*x.Tin2)...
            - ((v.F_outRP(:,i) + v.F_Rec(:,i)).*x.T2)...
            + (((p.UA_amb(i) - p.UA_RP(i))./p.C_p).*x.T2))./p.m_RPj...
            + x.Qrefr2(end) - x.Qamb2(end);

ddt.Qrefr2 = 0;
ddt.Qamb2  = 0;
ddt.Tin2   = 0;


A = [-(v.F_Rec(:,i)./p.m_RPj)-(v.F_outRP(:,i)./p.m_RPj)+(p.UA_amb(i)./(p.C_p*p.m_RPj))-(p.UA_RP(i)./(p.C_p*p.m_RPj)) 1 -1 ((v.F_outRP(:,i)+v.F_Rec(:,i))./p.m_RPj);    
       0 0 0 0;
       0 0 0 0;
       0 0 0 0]; 


ddt.P2 = A.*x.P2.*A' + p.Q2;
ddt.P2 = ddt.P2(1,1);

dxdt   = S2V(ddt, s.KFstatefields_RP(:,i));

