function o = TestIfOver(tStart, x, output, s, u, z, p, response, t,...
                        solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, plants, u_opt_CD)
    output.MV_Econ = plants;
    MV_SD = p.F_outSS_SD;
    MV_PT = p.F_outSS_PT;
    MV_CD = p.F_outSS_CD;
    MV1   = p.F_RecSS_RP;
    MV2   = p.F_RecSS_RP;
    MV3   = p.F_RecSS_RP;
    MV4   = p.F_RecSS_RP;
    MV5   = p.F_RecSS_RP;
    options = optimoptions('fmincon','Display','off');

 for i = ((tStart/p.Ts)+1):1:(tStart/p.Ts)+(p.Ts*p.N_Econ)
                v = RPIntermediatesKF(x, u, p, t(i), output, s, response);
                z = Meas(solSD, solPT, sol1, sol2, sol3, sol4, sol5, solCD, u, t(i), p, z, v);
                x = KalmanFilter(x, s, output, u, z, [t(i-1) t(i)], p, response);
                MV_PT(end+1) = 200 + 150*(plants(1)>=1) + 150*(plants(1)>=2) + 150*(plants(1)>=3) + 100*(plants(1)>=4) + 100*(plants(1)>=5);
                output.MV_PT = griddedInterpolant(t((tStart/p.Ts):i), MV_PT, 'previous');
                MV1(end+1) = MV1(1);
                output.MV1 = griddedInterpolant(t((tStart/p.Ts):i), MV1, 'previous');
                MV2(end+1) = MV2(1);
                output.MV2 = griddedInterpolant(t((tStart/p.Ts):i), MV2, 'previous');
                MV3(end+1) = MV3(1);
                output.MV3 = griddedInterpolant(t((tStart/p.Ts):i), MV3, 'previous');
                MV4(end+1) = MV4(1);
                output.MV4 = griddedInterpolant(t((tStart/p.Ts):i), MV4, 'previous');
                MV5(end+1) = MV5(1);
                output.MV5 = griddedInterpolant(t((tStart/p.Ts):i), MV5, 'previous');
                % Chill Dam
                u_opt_CD = fmincon(@(uMV) costCD(t(i), uMV, u, p, s, x, response, output), u_opt_CD, [], [], [], [],...
                           p.MV_min_CD, p.MV_max_CD, [], options);
                MV_CD(end+1) = u_opt_CD(1); 
                output.MV_CD = griddedInterpolant(t((tStart/p.Ts):i), MV_CD, 'previous');
                solCD = odextend(solCD, @(t,x) ChillDamODEs(s, p, x, u, t, output, response), t(i+1)); 
                response.CD(:,i) = deval(solCD, t(i));
                fprintf('i = %d\n',i); 
 end

 o = x.L_CD(end);