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

%v.T_inRP = ((v.T_CD.*output.MV_CD(t)) + (u.T_outPT(t).*u.F_outPT(t))) ./ (output.MV_CD(t) + u.F_outPT(t)); % We'd like to control this
                                                 