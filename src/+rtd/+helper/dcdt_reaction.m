function dcdt_reac = dcdt_reaction(c,k)

    dcdt_reac(1) = (- k(1)*c(1)*c(3) - k(4)*(c(6)+c(4))*c(1) + k(5)*c(6));
    dcdt_reac(2) = (-k(2)*c(2)*c(3) + k(3)*c(4)*c(5));
    dcdt_reac(3) = (- k(1)*c(1)*c(3) - k(2)*(c(7)+c(2))*c(3) + k(3)*(c(6)+c(4))*c(5));
    dcdt_reac(4) = (k(2)*c(3)*c(2)-k(3)*c(4)*c(5)-k(4)*c(4)*c(1));
    dcdt_reac(5) = (k(2)*c(3)*(c(2)+c(7)) - k(3)*(c(6)+c(4))*c(5));
    dcdt_reac(6) = (k(4)*c(4)*c(1)+k(2)*c(7)*c(3) - k(3)*c(6)*c(5));
    dcdt_reac(7) = (-k(2)*c(7)*c(3) + k(3)*c(6)*c(5) + k(1)*c(1)*c(3));
    dcdt_reac = dcdt_reac';

end