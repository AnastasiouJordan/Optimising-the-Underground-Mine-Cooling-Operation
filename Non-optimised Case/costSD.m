function J = costSD(tStart, u_MV, u, p, s, x, output)
        time = tStart + (0 : p.Ts : p.Ts*(p.N - 1)); % Prediction horizon length
        x0 = [x.L_SD(end); x.F_in_SD(end)]; % Set the starting point to the posterior state estimate from the KF

        output.MV_SD = griddedInterpolant(time, u_MV, 'previous');
  

        [~,x_response] = ode45(@(time,x) MPCSquareDamODEs(s,p,x,u,time,output), time, x0);
        x.L_SD = x_response(:,1)';
       
        p.SP_SD = p.L_SD_SP_SS; % Set the SP for the Square Dam Level to the output of the optimiser

        L_pred  = x.L_SD;     % Convert current inventory from mass to %
        SP_curr = p.SP_SD; % Convert current SP from mass to %

        % SET-POINT CONTROL
        %PV_cost = p.Q_Weight_SD*sum((SP_curr - L_pred).^2);
        MV_cost = p.R_Weight_SD*sum((u_MV(2:end) - u_MV(1:end-1)).^2);

        % LIMIT CONTROL
        PV_cost = p.Q_Weight_SD*sum(1./(L_pred - p.PV_min_SD).^4 + 1./(p.PV_max_SD - L_pred).^4);


 
        J = PV_cost + MV_cost;
end