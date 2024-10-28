function dxdt = KFSquareDamODEs(s, p, x, u, t, output)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
%x = V2S(x, s.KFstatefields);
x = V2S(x, s.KFstatefields_SD);


% Calculate state derivatives as structure
ddt.L_SD = (x.F_in_SD(end) - output.MV_SD(t) - p.m_evapSD)./p.m_SDmax*100; % Note: using u.control(1) is done
                                                           % to be able to
                                                           % specify it's
                                                           % values in the
                                                           % MAIN based on
                                                           % the control
                                                           % action
                                                           % calculated by
                                                           % the MPC.
                                                           % However it
                                                           % needs to be
                                                           % adjusted such
                                                           % that it can
                                                           % also be set to
                                                           % equal the
                                                           % original data
                                                           % set to
                                                           % calculate the
                                                           % original mass
                                                           % prediction.


ddt.F_in_SD = 0;

ddt.P_SD = p.A_SD.*x.P_SD.*p.A_SD' + p.Q_SD;
ddt.P_SD = ddt.P_SD(1,1);

dxdt = S2V(ddt,s.KFstatefields_SD);
