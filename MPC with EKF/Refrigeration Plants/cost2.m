function J = cost2(tStart, u_MV2, u, p, s, x, i)
        time = tStart + (0 : p.Ts : p.Ts*(p.N - 1)); % Prediction horizon length

        x0 = [x.T2(end); x.Qrefr2(end); x.Qamb2(end); x.Tin2(end)];
      
        output.MV2 = griddedInterpolant(time, u_MV2, 'previous');


        [~,x_response] = ode45(@(time,x) MPCFridgePlantsODEs2(p, x, u, time, output, s, i), time, x0);
        x.T2   = x_response(:,1)';
        x.Tin2 = x_response(:,4)';

        v.m_inRPtot = u.F_inRPtot(time); 
        v.n         = sum(u.s(time),2);
        u_s = u.s(time);
        if v.n > 0
            v.F_outRP = (u.F_outCD(time) + u.F_outPT(time))./(v.n.*(u_s(:,i)+0.001))';
        else
            v.F_outRP = u_s(:,i)+0.001;
        end
        v.F_Rec     = output.MV2(time);
        v.F_inRP    = v.F_outRP + v.F_Rec;       
        v.T_inRPtot = ((u.F_outCD(time).*u.T_outCD(time)) + ...
              (u.F_outPT(time).*u.T_outPT(time)))./...
              (u.F_outCD(time)+u.F_outPT(time)); 
        v.T_inRP    = ((v.F_Rec.*x.T2) + ...
              (v.F_outRP.*v.T_inRPtot))...
              ./ (v.F_Rec + v.F_outRP);


        T_pred = v.T_inRP;
        SP_curr = p.SP; 

        % Choose to track level SP or keep level within limits

%         % SP CONTROL
        %PV_cost = p.Q_Weight*sum(((SP_curr - L_pred).^2) + (12 - v.T_inRP).^2);
        PV_cost = u_s(1,i).*p.Q_Weight*sum((SP_curr - T_pred).^2);
        MV_cost = p.R_Weight*sum((u_MV2(2:end) - u_MV2(1:end-1)).^2);
        %LIMIT CONTROL
%         PV_cost = p.Q_Weight*sum(1./(L_pred - p.PV_min).^2 + 1./(p.PV_max - L_pred).^2);
%         MV_cost = p.R_Weight*sum((u_MV(2:end) - u_MV(1:end-1)).^2);

 
        J = PV_cost + MV_cost;
end