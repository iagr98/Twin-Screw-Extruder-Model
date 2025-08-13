function k = setup_kinparameters(T,ROH,kin)
% calculates the kinetic parameters of the reaction system for a given
% temperature T in C

    k_MC = rtd.helper.arrhenius(kin.kCM(1),kin.kCM(2),T);
    k_MC = k_MC + (kin.a) * ROH;
    % k_MC = 0;
    K_A = rtd.helper.arrhenius(kin.KA(1),kin.KA(2),T);
    k_A1 = 10e3;  
    % k_p = arrhenius(kin.kp(1),kin.kp(2),T);
    k_p = 300e3;
    M_eq = rtd.helper.calc_Meq(T);
    k_A2 = k_A1/K_A;
    k_d = k_p*M_eq;
    k_s = 1e6;
    k_te = rtd.helper.arrhenius(kin.kte(1),kin.kte(2),T);
    k_de = rtd.helper.arrhenius(kin.kde(1),kin.kde(2),T);
    k = [k_MC; k_A1; k_A2; k_p; k_d; k_s; k_te; k_de];

end
