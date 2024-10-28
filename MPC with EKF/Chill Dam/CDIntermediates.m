function v = CDIntermediates(x, u, p, t, output)
% Calculate intermediate process variables 
% (i.e. all variables which are neither exogeneous inputs 
%       nor state variables)
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs (measured variables)
%   p: structure of parameters
                                   
                                                        
v.F_RevFlow     = abs(min(u.F_outCD(t),0));               % Set the reversal of flow 
                                                          % to be the absolute value 
                                                          % of all the
                                                          % negative flows
                                                          % that leave the
                                                          % Chill Dam. 

% Calculate enthalpies, pressures & temperatures:
v.T_inCD = u.T_outRP(t);                                  % oC, Temperature into CD
v.T_CD   = (x.h_CD' - p.h_0 + (p.C_p * p.T_0)) ./ p.C_p;  % oC, Temperature out of CD, also measured
v.h_RP   = (p.C_p * (v.T_inCD - p.T_0)) + p.h_0;          % kJ/kg, Enthalpy of water into CD from the RPs
v.h_PT   = (p.C_p * (u.T_outPT(t) - p.T_0)) + p.h_0;      % kJ/kg, Enthalpy of water out of the PTs

v.T_inRP = ((v.T_CD.*output.MV(t)) + (u.T_outPT(t).*u.F_outPT(t))) ./ (output.MV(t) + u.F_outPT(t)); % We'd like to control this

% if v.m_2 > 0   
%     v.H_3 = x.h_CD.*v.m_3;
% else
%     v.H_3 = v.h_PT.*v.m_3;
% end
%v.H_3 = x.h_CD.*(v.m_3 > 0) + v.h_PT.*(v.m_3 < 0); % If the flowrate out of the CD is positive, let the enthalpy be equal to that of the water in the Chill Dam.
                                                   % If the flowrate out of the CD is negative, let the enthalpy be equal to that of the water leaving the Pre-Cooling Towers.