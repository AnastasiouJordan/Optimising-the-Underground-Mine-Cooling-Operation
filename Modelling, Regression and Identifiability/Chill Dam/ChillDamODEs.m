function dxdt = ChillDamODEs(s, p, x_vec, u, t)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x_vec, s.statefields);
v = CDIntermediates(x, u, p, t);


% if v.m_inCD > 3.760066666666660e+02
%     v.m_inCD = 3.760066666666660e+02;
% end
% if v.m_inCD < 36.087500000000254
%     v.m_inCD = 36.087500000000254;
% end
% 
% if v.m_outCD > 3.611966666666670e+02
%     v.m_outCD = 3.611966666666670e+02;
% end
% if v.m_outCD < 93.199999999999970
%     v.m_outCD = 93.199999999999970;
% end

% Calculate state derivatives as structure
%ddt.m_CD = u.F_out_PT_filtered(t) - u.F_UG_filtered(t) - u.F_Ice_filtered(t);
ddt.m_CD = v.m_1 - v.m_2 - v.m_3;

% ddt.h_CD = (v.H_1 - (v.m_1 .* x.h_CD) - ...
%             v.H_2 + (v.m_2 .* x.h_CD) - ...
%             v.H_3 + (v.m_3 .* x.h_CD) + ...
%             (p.UA_CD .* (u.T_Env(t) - v.T_CD))) ./ x.m_CD;

v.T_1 = v.T_inCD.*(v.m_1 > 0) + x.h_CD'.*(v.m_1 < 0); 
v.T_2 = x.h_CD'.*(v.m_2 > 0) + u.T_outPT(t).*(v.m_2 < 0); 
v.T_3 = x.h_CD';

ddt.h_CD = (v.T_1.*v.m_1 - ...
            v.T_2.*v.m_2 - ...
            v.T_3.*v.m_3 + ...
            (p.UA_CD .* (u.T_Env(t) - x.h_CD'))./p.C_p) ./ (x.m_CD);

% ddt.h_CD = (v.H_1 - (v.m_1 .* x.h_CD) - ...
%             v.H_2 + (v.m_2 .* x.h_CD) + ...  
%             (p.UA_CD .* (u.T_Env(t) - v.T_CD))) ./ x.m_CD;

% Map state derivative structure to vector
dxdt = S2V(ddt, s.statefields);
%dxdt = dxdt(:,1);