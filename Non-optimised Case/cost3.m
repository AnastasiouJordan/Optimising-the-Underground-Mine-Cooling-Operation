function J = cost3(tStart, u_MV3, u, p, s, x, i, output, response)
        time = tStart + (0 : p.Ts : p.Ts*(p.N - 1)); % Prediction horizon length

        x0 = [x.T3(end); x.Qrefr3(end); x.Qamb3(end); x.Tin3(end)];
      
        output.MV3 = griddedInterpolant(time, u_MV3, 'previous');

        v = RPIntermediatesKF(x, u, p, time, output, s, response);

        [~,x_response] = ode45(@(time,x) MPCFridgePlantsODEs3(p, x, u, time, output, s, i, v), time, x0);
        x.T3   = x_response(:,1)';
        x.Tin3 = x_response(:,4)';



        p.SP_RP = p.T_RP_SP_SS; 
        RPs_ON = output.MV_Econ(end);

        u_s = 1*(RPs_ON>2);

        T_pred = v.T_inRP(:,i);
        SP_curr = p.SP_RP;
        

%         % SP CONTROL
        %PV_cost = p.Q_Weight*sum(((SP_curr - L_pred).^2) + (12 - v.T_inRP).^2);
        PV_cost = u_s.*p.Q_Weight_RP*sum((SP_curr - T_pred).^2); 
        MV_cost = p.R_Weight_RP*sum((u_MV3(2:end) - u_MV3(1:end-1)).^2);
        %LIMIT CONTROL
%         PV_cost = p.Q_Weight_RP*sum(1./(L_pred - p.PV_min_RP).^2 + 1./(p.PV_max_RP - L_pred).^2);
%         MV_cost = p.R_Weight_RP*sum((u_MV(2:end) - u_MV(1:end-1)).^2);

 
        J = PV_cost + MV_cost;
end