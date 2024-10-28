function J = cost5(tStart, u_MV5, u, p, s, x, i, output, response)
        time = tStart + (0 : p.Ts : p.Ts*(p.N - 1)); % Prediction horizon length

        x0 = [x.T5(end); x.Qrefr5(end); x.Qamb5(end); x.Tin5(end)];
      
        output.MV5 = griddedInterpolant(time, u_MV5, 'previous');

        v = RPIntermediatesKF(x, u, p, time, output, s, response);

        [~,x_response] = ode45(@(time,x) MPCFridgePlantsODEs5(p, x, u, time, output, s, i, v), time, x0);
        x.T5   = x_response(:,1)';
        x.Tin5 = x_response(:,4)';


        p.SP_RP = p.T_RP_SP_SS; 
        RPs_ON = output.MV_Econ(end);

        u_s = 1*(RPs_ON>4);


        T_pred = v.T_inRP(:,i);
        SP_curr = p.SP_RP; % Convert current SP from mass to %


%         % SP CONTROL
        %PV_cost = p.Q_Weight*sum(((SP_curr - L_pred).^2) + (12 - v.T_inRP).^2);
        PV_cost = u_s.*p.Q_Weight_RP*sum((SP_curr - T_pred).^2); 
        MV_cost = p.R_Weight_RP*sum((u_MV5(2:end) - u_MV5(1:end-1)).^2);
        %LIMIT CONTROL
%         PV_cost = p.Q_Weight_RP*sum(1./(L_pred - p.PV_min_RP).^2 + 1./(p.PV_max_RP - L_pred).^2);
%         MV_cost = p.R_Weight_RP*sum((u_MV(2:end) - u_MV(1:end-1)).^2);

 
        J = PV_cost + MV_cost;
end