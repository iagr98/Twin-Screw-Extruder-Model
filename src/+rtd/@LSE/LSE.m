classdef LSE < rtd.RSE

    properties
    end

    methods
        function obj = LSE(extruder_geometry,screw,N)
            obj@rtd.RSE(extruder_geometry,screw,N);          
        end
        
     
        function obj = calcDFlow(obj)
            % Calculate the Drag Flow along the screw element accoring to
            % Booy and Vergnes
            
            %Calculation of current density and viscosity
            % obj = calcDensity(obj);
            obj.m_DFlow = obj.A1*obj.D_e^3*obj.N*obj.rho*obj.f;     %Caculation of Drag Flow
            obj.m_Flow(1,2) = obj.m_DFlow;
        end
        
        
    end
end