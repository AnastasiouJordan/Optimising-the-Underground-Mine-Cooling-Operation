function z = Meas(solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, u, t, p, z, v)


% Square Dam
z.L_SD(end+1) = solSD.y(end) + randn*sqrt(p.w_L_SD);
z.F_in_SD(end+1) = u.F_in_generatedSD(t) + randn*sqrt(p.w_F_in_SD);

% Pre-cooling Towers
z.L_PT(end+1) = solPT.y(1,end) + randn*sqrt(p.w_L_PT);
z.F_AntiSurge(end+1) = u.F_AntiSurge_generated(t) + randn*sqrt(p.w_F_Antisurge);
z.T_AntiSurge(end+1) = u.T_AntiSurge_generated(t) + randn*sqrt(p.w_T_Antisurge);
z.T_outSD(end+1) = u.T_out_generatedSD(t) + randn*sqrt(p.w_T_out_SD);
z.T_PT(end+1) = solPT.y(2,end) + randn*sqrt(p.w_T_PT);
z.T_amb(end+1) = u.T_amb_generated(t) + randn*sqrt(0.002);
z.H_Env(end+1) = u.H_Env_generated(t) + randn*sqrt(p.w_H_Env);


% Fridge Plants
z.T1(end+1) = sol1.y(end) + randn*sqrt(p.w_T_RP(1));
z.T2(end+1) = sol2.y(end) + randn*sqrt(p.w_T_RP(2));
z.T3(end+1) = sol3.y(end) + randn*sqrt(p.w_T_RP(3));
z.T4(end+1) = sol4.y(end) + randn*sqrt(p.w_T_RP(4));
z.T5(end+1) = sol5.y(end) + randn*sqrt(p.w_T_RP(5));
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
% z.Tin1(end+1) = v.T_inRP(end,1) + randn*sqrt(0.002);
% z.Tin2(end+1) = v.T_inRP(end,2) + randn*sqrt(p.w_T_in_RP(2));
% z.Tin3(end+1) = v.T_inRP(end,3) + randn*sqrt(p.w_T_in_RP(3));
% z.Tin4(end+1) = v.T_inRP(end,4) + randn*sqrt(p.w_T_in_RP(4));
% z.Tin5(end+1) = v.T_inRP(end,5) + randn*sqrt(p.w_T_in_RP(5));
z.Tin1(end+1) = v.T_inRP(end,1) + randn*sqrt(0.002);
z.Tin2(end+1) = v.T_inRP(end,2) + randn*sqrt(0.0001);
z.Tin3(end+1) = v.T_inRP(end,3) + randn*sqrt(0.0001);
z.Tin4(end+1) = v.T_inRP(end,4) + randn*sqrt(0.0001);
z.Tin5(end+1) = v.T_inRP(end,5) + randn*sqrt(0.0001);


% Chill Dam
z.L_CD(end+1) = solCD.y(1,end) + randn*sqrt(p.w_L_CD);
z.T_CD(end+1) = solCD.y(2,end) + randn*sqrt(p.w_T_CD);
z.F_Ice(end+1) = u.F_Ice_generated(t) + randn*sqrt(p.w_F_Ice);
z.F_UG(end+1) = u.F_UG_generated(t) + randn*sqrt(p.w_F_UG);