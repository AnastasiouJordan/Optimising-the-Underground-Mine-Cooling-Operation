function J = costCD(tStart, u_MV, u, p, s, x, response, output)
        time = tStart + (0 : p.Ts : p.Ts*(p.N - 1)); % Prediction horizon length
        x0 = [x.L_CD(end); x.T_CD(end); x.F_Ice(end); x.F_UG(end)]; % Set the starting point to the posterior state estimate from the KF

        output.MV_CD = griddedInterpolant(time, u_MV, 'previous');

        v = RPIntermediatesKF(x, u, p, time, output, s, response);
        v.T_PT = x.T_PT(end); v.T1 = x.T1(end); v.T2 = x.T2(end); v.T3 = x.T3(end); 
        v.T4 = x.T4(end); v.T5 = x.T5(end); v.T_amb = x.T_amb(end);
        
        [~,x_response] = ode45(@(time,x) MPCChillDamODEs(s,p,x,u,time,output,v), time, x0);
        x.L_CD = x_response(:,1)';
      
        p.SP_L_CD = p.L_CD_SP_SS; % Set the SP for the Chill Dam Level to the output of the optimiser
        p.SP_T_CD = p.T_CD_SP_SS; % Set the SP for the Chill Dam Outlet Temperature to the output of the optimiser

        L_pred    = x.L_CD;  % Convert current inventory from mass to %
        SP_curr_L = p.SP_L_CD;

        T_pred    = v.T_inRPtot*ones(1, length(time));  % Convert current inventory from mass to %
        SP_curr_T = p.SP_T_CD; % Convert current SP from mass to %

        PV_pred = [L_pred T_pred];
        SP_curr = [SP_curr_L SP_curr_T];
        % Choose to track level SP or keep level within limits

%         % SP CONTROL
        %PV_cost = 1.*sum(1./(T_pred - p.PV_min_T_CD).^2 + 1./(p.PV_max_T_CD - T_pred).^2) + 1.*sum(1./(L_pred - p.PV_min_L_CD).^4 + 1./(p.PV_max_L_CD - L_pred).^4);
        PV_cost = 0.01.*sum((SP_curr_T - T_pred).^2) + 1.*sum(1./(L_pred - p.PV_min_L_CD).^4 + 1./(p.PV_max_L_CD - L_pred).^4);
        MV_cost = p.R_Weight_CD.*sum((u_MV(2:end) - u_MV(1:end-1)).^2);
        %LIMIT CONTROL
       % PV_cost = p.Q_Weight_CD*sum(1./(L_pred - p.PV_min_CD).^4 + 1./(p.PV_max_CD - L_pred).^4);

 
        J = PV_cost + MV_cost;
end