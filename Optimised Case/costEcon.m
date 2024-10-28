function J = costEcon(tStart, MV_Econ, u, p, s, x, output)

        time = tStart + (0 : 3600 : 3600*(p.N_Econ - 1)); % Prediction horizon length


        %Econ_cost = 0.1.*sum(MV_Econ.*p.pull_full_efficiency.*p.whole_week_tariff_vec_minutes((time+28800)/60));
        Econ_cost = sum(MV_Econ.*p.pull_full_efficiency.*p.whole_week_tariff_vec_minutes(time/60));
        %MV_cost = sum((MV_Econ(2:end) - MV_Econ(1:end-1)).^2);
 
        J = Econ_cost;
end 