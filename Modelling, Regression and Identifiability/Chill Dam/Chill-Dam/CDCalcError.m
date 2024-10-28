function [E, x, v] = CDCalcError(pUAvec, u, p, s, t)

% Simulate the outputs of the system based on parameter estimates.
% First convert p_vector of inital parameter guesses to a structure:
p = V2S(pUAvec, p.regressedparameterfields, p);

x0.m_CD = u.L_CD(0)*p.m_CDmax/100;   
x0.h_CD = 20;          
x0_vec  = S2V(x0, s.statefields);

[~, x_vec] = ode45(@(t, x) ChillDamODEs(s, p, x, u, t), t, x0_vec);
x = V2S(x_vec', s.statefields);
v = CDIntermediates(x, u, p, t);

% Finally, we calculate the error based on this output of the simulation
% using the estimated parameters

%E = [v.T_CD - u.T_CD(t)];
%E = [v.L_CD - u.L_CD(t)];
E = [(v.T_CD - u.T_CD(t)) (v.L_CD - u.L_CD(t))];

