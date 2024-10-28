function J = costPT(tStart, u_MV, u, p, s, x, output)
        time = tStart + (0 : p.Ts : p.Ts*(p.N - 1)); % Prediction horizon length
        x0 = [x.L_PT(end); x.F_AntiSurge(end); x.T_AntiSurge(end); x.T_outSD(end); x.T_PT(end); x.T_amb(end); x.H_Env(end)]; % Set the starting point to the posterior state estimate from the KF

        output.MV_PT = griddedInterpolant(time, u_MV, 'previous');

        [~,x_response] = ode45(@(time,x) MPCPreCoolingTowerODEs(s,p,x,u,time,output), time, x0);
        x.L_PT = x_response(:,1)';
       

        p.SP_PT = p.L_PT_SP_SS; % Set the SP for the Pre-cooling Tower Level to the output of the optimiser

        L_pred  = x.L_PT;     % Convert current inventory from mass to %
        SP_curr = p.SP_PT; % Convert current SP from mass to %

        % Choose to track level SP or keep level within limits

%         % SP CONTROL
        %PV_cost = p.Q_Weight_PT*sum((SP_curr - L_pred).^2);
        MV_cost = p.R_Weight_PT*sum((u_MV(2:end) - u_MV(1:end-1)).^2);
        %LIMIT CONTROL
        PV_cost = p.Q_Weight_PT*sum(1./(L_pred - p.PV_min_PT).^4 + 1./(p.PV_max_PT - L_pred).^4);


 
        J = PV_cost + MV_cost;
end