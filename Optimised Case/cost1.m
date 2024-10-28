function J = cost1(tStart, u_MV1, u, p, s, x, i, output, response)
        time = tStart + (0 : p.Ts : p.Ts*(p.N - 1)); % Prediction horizon length

        x0 = [x.T1(end); x.Qrefr1(end); x.Qamb1(end); x.Tin1(end)];
     
        output.MV1 = griddedInterpolant(time, u_MV1, 'previous');

        v = RPIntermediatesKF(x, u, p, time, output, s, response);
        
        [~,x_response] = ode45(@(time,x) MPCFridgePlantsODEs1(p, x, u, time, output, s, i, v), time, x0);
        x.T1   = x_response(:,1)';
        x.Tin1 = x_response(:,4)';

     

        p.SP_RP = p.T_RP_SP_SS; 
        RPs_ON = output.MV_Econ(end);

        u_s = 1*(RPs_ON>=1);

        T_pred = v.T_inRP(:,i);
        SP_curr = p.SP_RP;
        
        
%         % SP CONTROL
        PV_cost = u_s.*p.Q_Weight_RP*sum((SP_curr - T_pred).^2);
        MV_cost = p.R_Weight_RP*sum((u_MV1(2:end) - u_MV1(1:end-1)).^2);
        %LIMIT CONTROL
%         PV_cost = p.Q_Weight_RP*sum(1./(L_pred - p.PV_min_RP).^2 + 1./(p.PV_max_RP - L_pred).^2);
%         MV_cost = p.R_Weight_RP*sum((u_MV(2:end) - u_MV(1:end-1)).^2);

 
        J = PV_cost + MV_cost;
end