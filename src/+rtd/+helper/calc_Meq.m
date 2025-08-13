function M_eq = calc_Meq(T)

    dHlc = -23.3e3;
    dSlc = -22;
%     dHlc = -50e3;
%     dSlc = -12;
    M_eq = exp(dHlc/(8.314*(T+273.15))-dSlc/8.314);

end