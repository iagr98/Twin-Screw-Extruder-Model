classdef Die < rtd.Reaktor
    %UNTITLED3 Summary of this class goes here
    %   Needs to be worked on. Initial assumption: Geometry is 2/3
    %   Kegelstumpf and 1/3 cylinder with die outlet diameter

    properties
        D_e       % outer diameter of screw
        D_i       % inner diameter of screw
        D_b       % barrel diameter
        C_l       % distance between screw centers
        d_in      % diameter of die at input
        d_out     % diameter of die at output
        d         % current diameter of die
        Die_count % Number of current reactor modelling the Die
        A_die     % coss-sectional area of die
        A_cr      % critical area = A_die for using same syntax of ohters reactors
        A_bm      % contact surface for heat transfer between barrel and melt
        k_p       % pressure loss factor               
        gamma_die % shear-thinning melt viscosity in the die                
        eta_av    % average viscosity of the melt in reaktor
        
        m_DFlow
        m_PFlow
        q_diss_ax
        r_max
        r_min
    end

    methods
        function obj = Die(extruder_geometry,screw,N, Die_count)  
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj@rtd.Reaktor(extruder_geometry,screw,N);            
            obj.D_e = extruder_geometry.D_e;
            obj.D_i = extruder_geometry.D_i;
            obj.D_b = extruder_geometry.D_b;
            obj.C_l = extruder_geometry.C_l;            
            obj.Die_count = Die_count;
            obj.L = screw.L/screw.n_reac;
            obj.d_in = screw.d_in;
            obj.d_out = screw.d_out;            
            obj.d = obj.d_in-0.5*(obj.d_in-obj.d_out)*obj.Die_count*obj.L/screw.L;
            obj.A_die = pi*obj.d^2/4;
            obj.r_max = 0.5*(obj.d + 0.5*(obj.d_in-obj.d_out)*obj.Die_count*obj.L/screw.L);
            obj.r_min = 0.5*(obj.d - 0.5*(obj.d_in-obj.d_out)*obj.Die_count*obj.L/screw.L);
            obj.A_cr = obj.A_die;
%             obj.A_bm = pi*obj.d*obj.L;  %Pending for correction
            obj.A_bm = pi*(obj.r_max+obj.r_min)*(obj.L*sqrt(1+((obj.r_max-obj.r_min)/screw.L)^2));
            obj.k_p = 128/pi*(obj.D_e/obj.d)^4;
            obj.V = obj.L*obj.d^2/4;
        end
        
        function obj = calcViscosityPLA(obj)
            % obj = calcViscosityPLA(obj) calculates viscosities of melt in
            % Die (eta_av) using Yasuda model(CY6)            
            if sum(obj.m) == 0
                obj.eta_av = obj.eta_min;  
                obj.gamma_die = 0;
            else
                if (obj.M_w~=0 && obj.w_p~=0 && obj.T_m(2)>obj.T_melt)
                    % m=0.395; eta_ref=3084; lambda_ref=0.0084; Ea=85300; T_ref=180+273; T_p=obj.T_m(2)+273;
                    % alfa=4.43; a=0.69; n=0.094; M_ref=100; R=8.314;

                    m=0.395; eta_ref=97; lambda_ref=0.0084; Ea=104e3; T_ref=180+273; T_p=obj.T_m(2)+273;
                    alfa=2.24; a=0.69; n=0.094; M_ref=100; R=8.314; n_w = 8;

                    volFlow = (sum(obj.m_Flow(:,1))+sum(obj.m_Flow(:,4))-sum(obj.m_Flow(:,3))-sum(obj.m_Flow(:,2)))/obj.rho;
                    volFlow(volFlow<0)=0;
                    obj.gamma_die = (3*(1-m)*8*volFlow)/((1-m)*obj.d^3*pi);
                    
                    obj.eta_av = eta_ref*power(obj.w_p,n_w)*(obj.M_w/M_ref)^alfa*exp((Ea/R)*(1/T_p-1/T_ref))*...
                        (1+(lambda_ref*(obj.M_w/M_ref)^alfa*exp((Ea/R)*(1/T_p-1/T_ref))*obj.gamma_die)^a)^((n-1)/a);


                    obj.eta_av(obj.eta_av<1e-3)=1e-3;
                else
                    obj.eta_av = obj.eta_min;
                end                    
            end            
        end

        function obj = calcDFlow(obj)
            % Needs a function to calc the drag flow. Drag flow in the die
            % is always 0
            obj.m_Flow(1,2:3) = 0;
        end

        function obj = calcPFlow(obj)            
            if obj.p(2) > obj.p(1) && obj.p(2) > obj.p(3)
                obj.m_Flow(2,2) = (obj.p(2)-obj.p(1))*obj.D_e^4*obj.rho/(obj.k_p*obj.eta_av*obj.L);
                obj.m_Flow(2,3) = (obj.p(2)-obj.p(3))*obj.D_e^4*obj.rho/(obj.k_p*obj.eta_av*obj.L);
            else
               if obj.p(2) > obj.p(3)
                   obj.m_Flow(2,3) = (obj.p(2)-obj.p(3))*obj.D_e^4*obj.rho/(obj.k_p*obj.eta_av*obj.L);
               end
               if obj.p(2) > obj.p(1)
                   obj.m_Flow(2,2) = (obj.p(2)-obj.p(1))*obj.D_e^4*obj.rho/(obj.k_p*obj.eta_av*obj.L);
               end
            end
            obj = calcViscosityPLA(obj);
        end   
        
        function obj = calcHeatFlow_bm(obj)
            % obj = calcHeatFlow_bm calculates heat flow between melt and barrel
            
            obj = calcViscosityPLA(obj);
            obj = calcHTC(obj);
            obj.HeatFlow_bm = obj.htc_b*obj.A_bm*(obj.T_b - obj.T_m(2));
        end
        
        function obj = calcheatFlow_diss(obj) 
            % obj = calcheatFlow_diss calculates heat flow by viscous dissipation per unit volume
            
            % ------------- Calculation of axial contribution -------------
            obj = calcViscosityPLA(obj);
            obj.q_diss_ax = (obj.D_e^4/(obj.k_p*obj.eta_av*obj.A_die))*((obj.p(3)-obj.p(1))/obj.L)^2;
            
            % Calculation of heat flow by viscous dissipation per unit volume
            obj.heatFlow_diss = obj.q_diss_ax*obj.m_Flow(1,3);
        end
    end
end