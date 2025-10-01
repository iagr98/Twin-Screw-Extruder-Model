classdef KB90 < rtd.KB45
    %Klasse zur Erstellung von Knetblcken

    properties
    end

    methods
        function obj = KB90(extruder_geometry,screw,N)
            obj@rtd.KB45(extruder_geometry,screw,N);
            obj.A0 = screw.A0;            
            obj.k_p = screw.A0;
        end

        function obj = calcDFlow(obj)
            %Calculation of current density and viscosity
            % obj = calcDensity(obj);
            obj.m_DFlow = 0;
            obj.m_Flow(1,3) = obj.m_DFlow;
        end
        
    end
end