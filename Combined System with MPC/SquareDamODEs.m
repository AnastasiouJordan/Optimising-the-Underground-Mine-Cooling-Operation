function dxdt = SquareDamODEs(s, p, x, u, t, output, z)
% Calculate the time-derivative of all state variables
%
% The function requires the following process variables as inputs:
%   t: time (scalar or vector)
%   x: structure of state variables
%   u: structure of exogeneous inputs
%   p: structure of parameters

% Map state vector to structure and calculate intermediate variables
x = V2S(x, s.statefields_SD);

% Calculate state derivatives as structure
ddt.L_SD = (u.F_in_generatedSD(t) - output.MV_SD(t) - p.m_evapSD)./p.m_SDmax*100; % Note: using u.control(1) is done
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

% Map state derivative structure to vector
dxdt = S2V(ddt,s.statefields_SD);
