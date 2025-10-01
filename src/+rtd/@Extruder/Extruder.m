classdef Extruder 
    properties
        extruder_geometry
        L_ges               % total length
        screw               % screw configuration
        AmountScrews        % number of different screw regions
        config
        reactorConfiguration
        volume
        m_Flow0             % Eingangsmassenstrom
        T                   % temperature
        P_U                 % ambient pressure
        n_discret = 0       % number of discretization points
        T_b
    end

    methods
        function obj = Extruder(D_e, D_i, D_b, C_l, T_b)
            obj.extruder_geometry.D_e = D_e; 
            obj.extruder_geometry.D_i = D_i; 
            obj.extruder_geometry.D_b = D_b; 
            obj.extruder_geometry.C_l = C_l; 
            load ('+rtd\+data\+variables\variables.mat', 'm_Flow0', 'P_U', 'T')
            obj.m_Flow0 = m_Flow0;
            obj.T = T;
            obj.P_U = P_U;
            obj.T_b = T_b;
        end
    
 
        function obj = designExtruder(obj,screw,N)
            
            % Calculate number of discrete points
            fields = fieldnames(screw);
            for i=1:length(fields)
                obj.n_discret = obj.n_discret + screw.(fields{i}).n_reac;
            end
            
            % Desing extruder
            obj.screw = struct2cell(screw);
            number = 1;
            Die_count = 1;
            obj.reactorConfiguration{1}.Die_count = Die_count;
            for ixScrew = 1:length(obj.screw)
                obj.volume = 0;   
                screw_i = obj.screw{ixScrew};
                for ixReactor = 1:screw_i.n_reac
                    switch screw_i.type
                        case 1
                            reactor = rtd.RSE(obj.extruder_geometry,screw_i,N);
                        case -1
                            reactor = rtd.LSE(obj.extruder_geometry,screw_i,N);
                        case 2                                                         
                            reactor = rtd.Die(obj.extruder_geometry,screw_i,N,Die_count);
                            Die_count = Die_count + 1; 
                        case 0
                            reactor = rtd.KB90(obj.extruder_geometry,screw_i,N);
                        case 3
                            reactor = rtd.KB45(obj.extruder_geometry,screw_i,N);
                    end
                    obj.reactorConfiguration{number}=reactor;
                    obj.reactorConfiguration{number}.number = number; 
                    obj.reactorConfiguration{number}.n_discret = obj.n_discret;                    
                    obj.reactorConfiguration{number}.T_b = obj.T_b(number);
                    obj.volume = obj.volume + obj.reactorConfiguration{number}.V; 
                    obj.reactorConfiguration{number}.f_ini = screw_i.f_ini;                    
                    number = number+1; 
                    reactor = [];
                end
            end
        end

        
        % function [obj,y0] = setInitials(obj)
        %     y0 = [];
        %     for i = 1:length(obj.reactorConfiguration)
        %         obj.reactorConfiguration{i} = obj.reactorConfiguration{i}.setInitialFilling;
        %         y_add = obj.reactorConfiguration{i}.m;
        %         y_add(8) = obj.T;
        %         y0 = [y0; y_add];
        %     end
        % end       
        
    end   
end