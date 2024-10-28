function z = Meas(sol, u, t, p, z)

% Please choose the measurements by commenting out
% the data you do not wish to use:

% Choose to take in Raw plant data:
% z.L(end+1)     = u.L_SD(t);
% z.F_in(end+1)  = u.F_inSD(t);


% Choose to take in Generated Noisy data:
z.L(end+1) = sol.y(1,end) + randn*sqrt(p.w_L);
z.F_in(end+1) = u.F_in_generated(t) + randn*sqrt(p.w_F_in);

