function c_start = setup_cstart(ratio)
% calculates initial molar concentrations c_start based on the ratio of 
% monomer : co-initiator : catalyst 
    
    M_LA = .14413;  % molecular weight of lactide in kg/mol
    rho = 1113.8;   % density of lactide at standard conditions
    M0 = rho/M_LA;  %concentration of Monomer in mol/m^3
    save variables/monomer0.mat M0
    % initial molar concentrations of reaction species in mol/m^3
    C0 = M0/ratio(1)*ratio(3);      % catalyst 
    ROH0 = M0/ratio(1)*ratio(2);    % Co-initiator
    A0 = 0;     % cation
    I0 = 0;
    D0 = 0;  % dormant chains
    R0 = 0;     % active chains
    
    c_start = [M0; ROH0; C0; I0; A0; R0; D0];
    
end