function z = Meas(sol1, sol2, sol3, sol4, sol5, u, t, p, z, v)

% Choose to take in Generated Noisy data:
z.T1(end+1) = sol1.y(end) + randn*sqrt(p.w_T(1));
z.T2(end+1) = sol2.y(end) + randn*sqrt(p.w_T(2));
z.T3(end+1) = sol3.y(end) + randn*sqrt(p.w_T(3));
z.T4(end+1) = sol4.y(end) + randn*sqrt(p.w_T(4));
z.T5(end+1) = sol5.y(end) + randn*sqrt(p.w_T(5));

v.Q_refr  = u.Q_refr_generated(t);
v.Q_amb   = u.Q_amb_generated(t);


z.Qrefr1(end+1) = v.Q_refr(:,1) + randn*sqrt(p.w_Q_refr(1));
z.Qrefr2(end+1) = v.Q_refr(:,2) + randn*sqrt(p.w_Q_refr(2));
z.Qrefr3(end+1) = v.Q_refr(:,3) + randn*sqrt(p.w_Q_refr(3));
z.Qrefr4(end+1) = v.Q_refr(:,4) + randn*sqrt(p.w_Q_refr(4));
z.Qrefr5(end+1) = v.Q_refr(:,5) + randn*sqrt(p.w_Q_refr(5));

p.w_Q_amb = 5E-9;
z.Qamb1(end+1) = v.Q_amb(:,1) + randn*sqrt(p.w_Q_amb);
z.Qamb2(end+1) = v.Q_amb(:,2) + randn*sqrt(p.w_Q_amb);
z.Qamb3(end+1) = v.Q_amb(:,3) + randn*sqrt(p.w_Q_amb);
z.Qamb4(end+1) = v.Q_amb(:,4) + randn*sqrt(p.w_Q_amb);
z.Qamb5(end+1) = v.Q_amb(:,5) + randn*sqrt(p.w_Q_amb);

z.Tin1(end+1) = v.T_inRP(end,1) + randn*sqrt(0.002);
z.Tin2(end+1) = v.T_inRP(end,2) + randn*sqrt(p.w_T_in(2));
z.Tin3(end+1) = v.T_inRP(end,3) + randn*sqrt(p.w_T_in(3));
z.Tin4(end+1) = v.T_inRP(end,4) + randn*sqrt(p.w_T_in(4));
z.Tin5(end+1) = v.T_inRP(end,5) + randn*sqrt(p.w_T_in(5));