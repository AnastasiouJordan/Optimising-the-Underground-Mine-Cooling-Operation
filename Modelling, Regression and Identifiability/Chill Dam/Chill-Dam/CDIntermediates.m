function v = CDIntermediates(x, u, p, t)
% Calculate intermediate process variables 
% (i.e. all variables which are neither exogeneous inputs 
%       nor state variables)
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs (measured variables)
%   p: structure of parameters


% Calculate mass & volumetric flowrates:
v.m_1       = u.F_in_filtered(t).*p.rho_Water./1000;  % kg/s, Inlet mass flowrate based on filtered measurements
v.m_2       = u.F_out_filtered(t).*p.rho_Water./1000; % kg/s, Outlet mass flowrate based on filtered measurements
v.m_3       = u.F_Ice_filtered(t).*p.rho_Water./1000; % kg/s, Outlet Ice plant flowrate based on filtered steps                                                     
                                                        
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

v.L_CD   = x.m_CD' ./ p.m_CDmax * 100;                    % %, Level prediction in CD based on the mass balance ODE

% Determine the enthalpies that will be assigned to the inlet and outlet
% streams
v.H_1 = v.h_RP.*(v.m_1 > 0) + x.h_CD'.*(v.m_1 < 0); % If the flowrate into the CD is positive, let the enthalpy be equal to that of the water leaving the Fridge Plants.
                                                   % If the flowrate into the CD is negative, let the enthalpy be equal to that of the water in the Chill Dam.
v.H_2 = x.h_CD'.*(v.m_2 > 0) + v.h_PT.*(v.m_2 < 0); % If the flowrate out of the CD is positive, let the enthalpy be equal to that of the water in the Chill Dam.
                                                   % If the flowrate out of the CD is negative, let the enthalpy be equal to that of the water leaving the Pre-Cooling Towers.
v.H_3 = x.h_CD'.*v.m_3;
% if v.m_2 > 0   
%     v.H_3 = x.h_CD.*v.m_3;
% else
%     v.H_3 = v.h_PT.*v.m_3;
% end
%v.H_3 = x.h_CD.*(v.m_3 > 0) + v.h_PT.*(v.m_3 < 0); % If the flowrate out of the CD is positive, let the enthalpy be equal to that of the water in the Chill Dam.
                                                   % If the flowrate out of the CD is negative, let the enthalpy be equal to that of the water leaving the Pre-Cooling Towers.