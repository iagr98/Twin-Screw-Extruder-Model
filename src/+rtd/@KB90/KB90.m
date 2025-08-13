classdef KB90 < rtd.Reaktor
    %Klasse zur Erstellung von Knetblcken

    properties
        psi       % 
        alpha     % tip angle
        D_e       % outer diameter of screw
        D_i       % inner diameter of screw
        D_b       % barrel diameter
        C_l       % distance between screw centers
        theta     % staggering angle of discs
        n_D       % number of discs
        t_eq      % equivalent pitch 
        n_F       % number of flights
        h         % gap height
        H         % channel height
        H_av      % average shannel height
        L_KB      % total length of KB section
        A_Scr     % cross-sectional area of screw
        L_Scr     % circumference of screw
        A_cr      % cross-sectional area of twin screw
        A_Surf    % surface area of conveying element
        A_bm      % contact surface for heat transfer between barrel and melt
        A0        % parameter for pressure flow calculation
        A1        % parameter for flow calculation
        A2        % parameter for flow calculation
        A3        % parameter for flow calculation
        L_D       % length of a disk (70% of L_Disk is L_D)
        L_off     % length of offset (30% of L_Disk is L_D)
        eps_ga    % volume fraction of the gap
        eps_off   % volume fraction of the offset
        eps_ch    % volume fraction of the channel
        gamma_ga  % shear-thinning melt viscosity in the gap
        gamma_ch  % shear-thinning melt viscosity in the channel
        gamma_off % shear-thinning melt viscosity in the offset between disks        
        eta_ga    % viscosity of the melt in the gap
        eta_off   % viscosity of the melt in the offset
        eta_ch    % viscosity of the melt in the channel
        eta_av    % average viscosity of the melt in reaktor
        
        k_d = 0
        k_p = 0
        m_DFlow =0
        m_PFlow = 0 
        q_diss_ax = 0
        q_diss_circ = 0
    end

    methods
        function obj = KB90(extruder_geometry,screw,N)
            %KB Construct an instance of class kneading block

            obj@rtd.Reaktor(extruder_geometry,screw,N);

            obj.D_e = extruder_geometry.D_e;
            obj.D_i = extruder_geometry.D_i;
            obj.D_b = extruder_geometry.D_b;
            obj.C_l = extruder_geometry.C_l;
            obj.theta = screw.theta;
            obj.n_D = screw.n_D;
            obj.t_eq = [];
            obj.n_F = screw.n_F;
            obj.L = screw.L/screw.n_reac;
            obj.L_KB = screw.L;
            obj.h = (obj.D_b - obj.D_e)/2;
            obj.H = (obj.D_b - obj.D_i)/2;
            obj.L_D = 0.7*(obj.L_KB/obj.n_D);
            obj.L_off = 0.3*(obj.L_KB/obj.n_D);

            obj.psi = acos(obj.C_l/obj.D_e);
            obj.alpha = pi/obj.n_F-2*obj.psi;
            obj.A_Scr = ((obj.D_i+obj.D_e)/2*((obj.D_i+obj.D_e)/2*obj.psi-...
                obj.D_e/2*sin(obj.psi))+obj.alpha/8*(obj.D_e^2+obj.D_i^2))*obj.n_F;
            obj.L_Scr = (obj.D_e+obj.D_i)/2*pi;
            obj.A_cr = 1/2*((pi-obj.psi)*obj.D_b^2+obj.C_l*sqrt(obj.D_b^2-obj.C_l^2))-2*obj.A_Scr;
            obj.A_bm = 2*obj.D_b*(pi-obj.psi)*obj.L;
            
            obj.eps_ga = obj.alpha*obj.D_e*obj.n_F*obj.h*obj.L_D*obj.n_D/(obj.A_cr*obj.L_KB);
            obj.eps_off = ((obj.alpha*obj.D_e*(obj.D_e-obj.C_l)*obj.L_off*(obj.n_D-1))/(2*obj.A_cr*obj.L_KB))...
                *(obj.psi*obj.n_F/(2*pi))*abs(obj.theta/90);
            obj.eps_ch = 1 - obj.eps_ga - obj.eps_off;
            obj.H_av = (obj.D_b^2*0.25*pi-obj.A_Scr-obj.alpha*0.5*obj.D_e*obj.h*obj.n_F)/...
                (0.5*obj.D_b*(4*obj.psi+obj.alpha)*obj.n_F);
            obj.gamma_ga = pi*obj.D_e*obj.N/obj.h;
            obj.gamma_ch = pi*obj.D_i*obj.N/obj.H_av;
            obj.gamma_off = 2*obj.C_l*pi*obj.N/obj.L_off;

            obj.A0 = screw.A0;            
            obj.k_p = screw.A0;
            obj.V = obj.A_cr*obj.L;
                      

        end
        
        function obj = calcViscosityPLA(obj)
            % obj = calcViscosityPLA(obj) calculates viscosities of melt in
            % gap, channel and offset (eta_ga, eta_ch, eta_off) and the 
            % average viscosity in (eta_av)in the reactor using Yasuda model(CY6)            
            if sum(obj.m) == 0
                obj.eta_av = obj.eta_min; obj.eta_ga = obj.eta_min; obj.eta_ch = obj.eta_min; obj.eta_off = obj.eta_min;  % update             
            else 
                if (obj.M_w~=0 && obj.w_p~=0 && obj.T_m(2)>obj.T_melt)
                    % eta_ref=3084; lambda_ref=0.0084; Ea=85300; T_ref=180+273; T_p=obj.T_m(2)+273;
                    % alfa=4.43; a=0.69; n=0.094; M_ref=100; R=8.314; 

                    eta_ref=97; lambda_ref=0.0084; Ea=104e3; T_ref=180+273; T_p=obj.T_m(2)+273;
                    alfa=2.24; a=0.69; n=0.094; M_ref=100; R=8.314; n_w = 8;
                    
                    obj.eta_ga = eta_ref*power(obj.w_p,n_w)*(obj.M_w/M_ref)^alfa*exp((Ea/R)*(1/T_p-1/T_ref))*...
                        (1+(lambda_ref*(obj.M_w/M_ref)^alfa*exp((Ea/R)*(1/T_p-1/T_ref))*obj.gamma_ga)^a)^((n-1)/a);
                    obj.eta_ch = eta_ref*power(obj.w_p,n_w)*(obj.M_w/M_ref)^alfa*exp((Ea/R)*(1/T_p-1/T_ref))*...
                        (1+(lambda_ref*(obj.M_w/M_ref)^alfa*exp((Ea/R)*(1/T_p-1/T_ref))*obj.gamma_ch)^a)^((n-1)/a);               
                    obj.eta_off = eta_ref*power(obj.w_p,n_w)*(obj.M_w/M_ref)^alfa*exp((Ea/R)*(1/T_p-1/T_ref))*...
                        (1+(lambda_ref*(obj.M_w/M_ref)^alfa*exp((Ea/R)*(1/T_p-1/T_ref))*obj.gamma_off)^a)^((n-1)/a);               
                    obj.eta_av = (obj.eta_ga*obj.eps_ga + obj.eta_off*obj.eps_off + obj.eta_ch*obj.eps_ch*obj.f)/...
                        (obj.eps_ga + obj.eps_off + obj.eps_ch*obj.f);

                            
                    obj.eta_av(obj.eta_av<1e-3)=1e-3;
                    obj.eta_ga(obj.eta_ga<1e-3)=1e-3;
                    obj.eta_ch(obj.eta_ch<1e-3)=1e-3;
                    obj.eta_off(obj.eta_off<1e-3)=1e-3;
                else
                    obj.eta_av = obj.eta_min;
                    obj.eta_ga = obj.eta_min;
                    obj.eta_ch = obj.eta_min;
                    obj.eta_off = obj.eta_min;
                end
            end            
        end

        function obj = calcDFlow(obj)
            
            %Calculation of current density and viscosity
            % obj = calcDensity(obj);
            obj.m_DFlow = 0;
            obj.m_Flow(1,3) = obj.m_DFlow;

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
            obj.q_diss_ax = (obj.D_e^4/(obj.k_p*obj.eta_av*obj.A_cr))*((obj.p(3)-obj.p(1))/obj.L)^2;
            
            % -------- Calculation of circumferential contribution --------           
            q_diss_ga = obj.eta_ga*obj.gamma_ga^2;
            q_diss_off = obj.eta_off*obj.gamma_off^2;
            q_diss_ch = 4*obj.eta_ch*obj.gamma_ch^2;
            
            obj.q_diss_circ = (q_diss_ga*obj.eps_ga + q_diss_off*obj.eps_off + q_diss_ch*obj.eps_ch*obj.f)/...
                (obj.eps_ga + obj.eps_off + obj.eps_ch*obj.f);
           
            
            % Calculation of heat flow by viscous dissipation per unit volume
            obj.heatFlow_diss = obj.q_diss_circ + obj.q_diss_ax*obj.m_Flow(1,3);
        end
    end
end