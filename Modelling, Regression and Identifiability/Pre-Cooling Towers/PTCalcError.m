function [E, x, v] = PTCalcError(pmAirVec, u, p, s, t)

% Simulate the outputs of the system based on parameter estimates.
% First convert p_vector of inital parameter guesses to a structure:

 
p = V2S(pmAirVec, p.regressedparameterfields, p);

x0.m_PT = u.L_PTA(0).*(p.m_PTAmax + p.m_PTBmax)/100;
x0.h_PT = 80;                     % kJ/kgK, Initial value for enthalpy
                                  % of water leaving PTs
x0_vec  = S2V(x0, s.statefields);

[~, x_vec] = ode45(@(t, x) PreCoolingTowersODEs(s, p, x, u, t), t, x0_vec);
x = V2S(x_vec', s.statefields);
v = PTIntermediates(x, u, p, t);

% Finally, we calculate the error based on this output of the simulation
% using the estimated parameters

E = [(v.T_PT - u.T_PT(t)) (v.L_PT - u.L_PTA(t))];
