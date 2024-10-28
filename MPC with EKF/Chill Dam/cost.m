function J = cost(tStart, u_MV, u, p, s, x)
        time = tStart + (0 : p.Ts : p.Ts*(p.N - 1)); % Prediction horizon length
        x0 = [x.L(end); x.F_in(end); x.h_CD(end)]; % Set the starting point to the posterior state estimate from the KF

        output.MV = griddedInterpolant(time, u_MV, 'previous');


        % Set the inlet flow to remain constant over the prediction horizon
%         x.F_in = x.F_in(end)*ones(1,p.N); % Creates a vector of the current inlet flowrate of length N
%         x.F_in = griddedInterpolant(time, x.F_in); % Converts it to the format required by the ODE: a gridded
%                                                        % interpolant of the
%                                                        % same inlet
%                                                        % flowrate over the
%                                                        % entire prediction
%                                                        % horizon length

        [~,x_response] = ode45(@(time,x) MPCChillDamODEs(s,p,x,u,time,output), time, x0);
        x.L = x_response(:,1)';
       



        L_pred  = x.L;     % Convert current inventory from mass to %
        SP_curr = p.SP; % Convert current SP from mass to %

        % Choose to track level SP or keep level within limits

%         % SP CONTROL
        %PV_cost = p.Q_Weight*sum(((SP_curr - L_pred).^2) + (12 - v.T_inRP).^2);
        PV_cost = p.Q_Weight*sum((SP_curr - L_pred).^2); 
        MV_cost = p.R_Weight*sum((u_MV(2:end) - u_MV(1:end-1)).^2);
        %LIMIT CONTROL
%         PV_cost = p.Q_Weight*sum(1./(L_pred - p.PV_min).^2 + 1./(p.PV_max - L_pred).^2);
%         MV_cost = p.R_Weight*sum((u_MV(2:end) - u_MV(1:end-1)).^2);

 
        J = PV_cost + MV_cost;
end