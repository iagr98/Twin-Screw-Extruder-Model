function [T_b, N_tot] = calcBarrelTemperature(screw)
    screw_cell = struct2cell(screw);
    N_tot = 0;                         % total number of reactors
    idx = 0;
    for i = 1:length(screw_cell)
        N_tot = N_tot + screw_cell{i}.n_reac; 
        n = screw_cell{i}.n_reac;
        T_b(idx+1 : idx+n) = screw_cell{i}.T_barrel*ones(1,n);
        idx = idx + n;
    end
end