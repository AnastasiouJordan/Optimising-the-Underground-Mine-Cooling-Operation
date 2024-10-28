function [E, x, v] = SDCalcError(pmEvapVec, u, p, s, t)

% Simulate the outputs of the system based on parameter estimates.
% First convert p_vector of inital parameter guesses to a structure:
p = V2S(pmEvapVec, p.regressedparameterfields, p);


x0.m_SD = u.L_SD(0)*p.m_SDmax/100;         
x0_vec  = S2V(x0, s.statefields);

[~, x_vec] = ode45(@(t, x) SquareDamODEs(s, p, x, u, t), t, x0_vec);
x = V2S(x_vec', s.statefields);
v = SDIntermediates(x, u, p, t);

% Finally, we calculate the error based on this output of the simulation
% using the estimated parameters

E = (v.L_SD - u.L_SD(t));
