function k = arrhenius(Ea,k0,T)

    R = 8.314;
    k = k0*exp(-Ea/(R*T));

end