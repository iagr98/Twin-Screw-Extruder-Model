classdef ScrewElement < matlab.mixin.Heterogeneous

    properties
        % Prozessbezogene Gren
        m_Flow0
        m_Flow = [0 0 0 0; 0 0 0 0] %Flow-Matrix im Schema [Dragflow in from left,
        % Dragflow out to left, Dragflow out to right, Dragflow in from
        % right; Pressureflow in from left, Pressureflow out to left,
        % Pressureflow out to right,
        % Pressureflow in from right]
        p = [1e5 1e5 1e5]   %Druck-Matrix im Schema [Reaktor i-1, Reaktor i, Reaktor i+1]
        f = 0               %Fllstand definiert durch die Kanalbreite
        f_ini               %Anfangsfllstand
        T = 0               %Temperatur im Reaktor
        T_b                 %Temperature of barrel
        T_m = [0 0 0]       %Temperature of the melt [T_m(i-1) T_m(i) T_m(i+1)]
        T_melt              %Melting temperature of polymer [C]
        gamma = 0
        lambda = 0.2        %Heat conductivity
        screw
        w = zeros(7,3)      %[Monomer, ROH, Catalyst, I, H+, R, D]
        rho
        eta   
        eta_new
        c_p                 % heat capacity of polymer [J/kgK]
        mol = zeros(7,1)    %[Monomer, Co-Initiator, Catalyst, Initiator, Acid, R, D] / mol
        M0 = 7.7277e3;      % initial monomer concentration / mol/m3
        M_rep = 72e-3;      % molar mass of the monomeric repeating unit / kg/mol
        h_R = -20e3         % Enthalpy of reaction / J/mol
        time = 0
        m = 0
        conv =0
        molarMass = zeros(7,1)
        c = zeros(7,3)      % Concentration of substance / mol/L
        Mn = zeros(3,1)
        dydt                % ode equations for current reactor
        N                   % Drehzal
        k                   % kinetic parameters of reaction
        ratio
        eta_params          % parameter for calculation of the viscosity
        is_reaction         % can be 'on' or 'off'
        number
        n_discret
        L                   % length of reactor
        V                   % maximum volume of reactor
        htc_b               % Heat transfer coefficient between barrel-melt
        HeatFlow_bm         % Heat flow between barrel and melt
        heatFlow_diss       % Heat flow per unit volume due to viscous dissipation
        melt_length         % Length of melting area in m
        
        term0 = 0
        term1 = 0
        term2 = 0
        term3 = 0
        term4 = 0
        dTdt = 0

        dmdt_flow
        dmdt_reac
        dcdt_reac = zeros(7,1)
        C1_kin
        ROH0
        dmdt_analysis = zeros(7,1)
        dVdt = 0        
        M_w = 10
        w_p = 1
        eta_min = 0.1
        flag_term4 = 0

                
    end
    methods
        %Konstruktor der Klasse
        function obj = ScrewElement(extruder_geometry,screw,N)            

            obj.f = screw.f_ini;
            % obj.T_b = screw.T_b;
            obj.N = N;
            obj.rho = 1120;
            obj.c_p = 2260;
            
            load +rtd\+data\+variables\variables.mat T m_Flow0 T_melt
            load +rtd\+data\+variables\kinetic.mat M C1_kin ratio
            load +rtd\+data\+variables\viscosity.mat

            obj.eta_params.B = B;
            obj.eta_params.E_A = E_a;
            obj.eta_params.R = R;
            obj.eta_params.a = a;
            obj.eta_params.eta_0_ref = eta_0_ref;

            obj.ratio = ratio;

            obj.ROH0 = obj.M0/obj.ratio(1)*obj.ratio(2);
            
            obj.m_Flow0 = m_Flow0;
            
            obj.molarMass = M;

            obj.T = T;
            obj.T_melt = T_melt;
            % obj = calcDensity(obj);
            obj.screw = screw;   
            obj.C1_kin = C1_kin;
            obj.k = rtd.Reactor.setup_kinparameters(150,obj.ROH0,obj.C1_kin);
            obj.eta = 100;
            obj.melt_length = 0;

            volume = sum(obj.m)/obj.rho;
            obj.conv = (obj.M0*volume-obj.mol(1))/(obj.M0*volume);
            obj.conv(obj.conv<0)=0;
                    
        end

        function obj = calcDensity(obj)  
            % hier muss noch genauer werden, Stand jetzt wird nur der
            % Wert fr 150C bei den properties vorgegeben
            rho_PLA = 1.1452/(1+0.00074*((obj.T-273.15)-150));
            obj.rho = rho_PLA*1e3;            
        end
        
        function obj = calcHTC(obj)
            % obj = calcHTC(obj) calculates heat transfer coefficient
            % between melt and barrer (htc_b)
            
            if obj.screw.type == 2
                D = obj.d;          %Pending for right definition
            else
                D = obj.D_b;
            end
            obj = calcViscosityPLA(obj);
            Re = obj.N*D^2*obj.rho/obj.eta_av;
            Pr = 500;
            Nu = 26.8*Re^0.838*Pr^0.33;
            obj.lambda = 0.2;
            obj.htc_b = Nu*obj.lambda/D;
        end

        function obj = calcTerm4(obj)
            mass_flow_in = sum(obj.m_Flow(:,1)) + sum(obj.m_Flow(:,4));
            tau = mass_flow_in /sum(obj.m);
            dc_dt = (obj.c(1,3)-obj.c(1,2))/tau;
            if (obj.number==obj.n_discret)
                dc_dt = 0;
            end
            obj.term4 = -dc_dt*obj.h_R/(obj.f*obj.rho*obj.c_p);
        end
        
        function obj = ode(obj)
            % ode sets up differntial equations for CSTR: dm/dt & dT/dt
            % obj = ode(obj) returns the DEs describing reactors behavior
            
            % ------------------------- MASS BALANCES -------------------------
            % dm/dt = dm/dt_flow + dm/dt_reac
            % [M,ROH,C,I,H+,R,D]    

           obj.dmdt_flow = sum(obj.m_Flow(:,1))-sum(obj.m_Flow(:,2))...
               -sum(obj.m_Flow(:,3))+sum(obj.m_Flow(:,4));
           obj.dmdt_flow(isnan(obj.dmdt_flow)) = 0;
           obj.dmdt_flow(obj.dmdt_flow<0)=0;
           obj.dydt(1) = obj.dmdt_flow;           
           
           % ------------------------- ENERGY BALANCES -------------------------
           % dT/dt = term0 + term1 + term2 - term3
           
           f_temp = max(obj.f,0.01);
           if (sum(obj.m) > 1e-10)
               obj.term0 =(sum(obj.m_Flow(:,1))*obj.T_m(1)*obj.c_p + sum(obj.m_Flow(:,4))*obj.T_m(3)*obj.c_p - sum(sum(obj.m_Flow(:,2:3)))*obj.T_m(2)*obj.c_p )...
                   /(f_temp*obj.rho*obj.V*obj.c_p);
               obj.term1 = (obj.f == 0) * 0 + (obj.f ~= 0) * obj.HeatFlow_bm/(f_temp*obj.rho*obj.V*obj.c_p);
               obj.term2 = (obj.f == 0) * 0 + (obj.f ~= 0) *obj.heatFlow_diss/(obj.rho*obj.c_p);
               obj.term3 = obj.T_m(2)*(sum(obj.m_Flow(:,1)) + sum(obj.m_Flow(:,4)) - sum(obj.m_Flow(:,3)) - sum(obj.m_Flow(:,2)))/(f_temp*obj.rho*obj.V); 
               if (obj.flag_term4 == 1), obj = calcTerm4(obj); end
    
               obj.dTdt = obj.term0 + obj.term1 + obj.term2 - obj.term3 + obj.term4;
               
           else
               obj.dTdt = 0;
           end
           obj.dTdt(isnan(obj.dTdt)|isinf(obj.dTdt)) = 0;           
           obj.dTdt(obj.dTdt<0)=0;
           obj.dydt(2) = obj.dTdt;
           
           % -------------------------- REACTION ON --------------------------
           
           % if strcmp(obj.is_reaction,'on') && (obj.T_m(2)>obj.T_melt)
           %     obj.mol = obj.m./obj.molarMass;
           %     obj.mol(isnan(obj.mol)) = 0;
           %     volume = obj.V*obj.f*1000; % kinetic parameters are in liter
           %     c = obj.mol./volume;
           %     c(isnan(c) | isinf(c)) = 0;
           %     % [k_MC,k_A1,k_A2,k_p,k_d] = deal(obj.k(1),obj.k(2),obj.k(3),obj.k(4),obj.k(5));
           %     % [M,ROH,C,I,A,R,D] = deal(c(1),c(2),c(3),c(4),c(5),c(6),c(7));
           %     obj.dcdt_reac(1) = (- obj.k(1)*c(1)*c(3) - obj.k(4)*(c(6)+c(4))*c(1) + obj.k(5)*c(6));
           %     obj.dcdt_reac(2) = (-obj.k(2)*c(2)*c(3) + obj.k(3)*c(4)*c(5));
           %     obj.dcdt_reac(3) = (- obj.k(1)*c(1)*c(3) - obj.k(2)*(c(7)+c(2))*c(3) + obj.k(3)*(c(6)+c(4))*c(5));
           %     obj.dcdt_reac(4) = (obj.k(2)*c(3)*c(2)-obj.k(3)*c(4)*c(5)-obj.k(4)*c(4)*c(1));
           %     obj.dcdt_reac(5) = (obj.k(2)*c(3)*(c(2)+c(7)) - obj.k(3)*(c(6)+c(4))*c(5));
           %     obj.dcdt_reac(6) = (obj.k(4)*c(4)*c(1)+obj.k(2)*c(7)*c(3) - obj.k(3)*c(6)*c(5));
           %     obj.dcdt_reac(7) = (-obj.k(2)*c(7)*c(3) + obj.k(3)*c(6)*c(5) + obj.k(1)*c(1)*c(3));
           %     % dVdt = sum(obj.dmdt_flow)/obj.rho*1000; % in L/s
           %     obj.dVdt = 0;
           %     obj.dmdt_reac = obj.dcdt_reac.*obj.molarMass*volume+obj.dVdt*c.*obj.molarMass;
           %     obj.dydt(1:7) = obj.dmdt_flow + obj.dmdt_reac;
           %     dMdt = -obj.M_rep/(obj.M0/500)*obj.dcdt_reac(1)*1000; % in kg/mol, Nherung fr 
           %     % instantane Aktivierung und kein Kettenstart durch Kat,
           %     % noch allgemeiner aufschreiben, muss evtl mit in ode
           %     % dMdt = 2*obj.M_rep*((dcdt_reac(7)*(M-obj.M0/1000)-D*dcdt_reac(1)...
           %     %     -R*dcdt_reac(1)+M*dcdt_reac(6)-obj.M0/1000*dcdt_reac(6))/((D+R)^2));
           %     dMdt(isnan(dMdt)||isinf(dMdt)) = 0;
           %     obj.dydt(6) = obj.dydt(6) + dMdt*c(6)*volume;
           %     obj.dydt(7) = obj.dydt(7) + dMdt*c(7)*volume;
           %     obj.dmdt_analysis = obj.dmdt_reac;
           %     obj.dmdt_analysis(6) = obj.dmdt_reac(6) + dMdt*c(6)*volume;
           %     obj.dmdt_analysis(7) = obj.dmdt_reac(7) + dMdt*c(7)*volume;
           %     obj.term4 = obj.dcdt_reac(1)*obj.h_R/(f_temp*obj.rho*obj.c_p);
           %     obj.dydt(8) = obj.dydt(8) + obj.term4;
           %     obj.dydt(isnan(obj.dydt)) = 0; 
           % end
           
           
        end
        
        function obj = calcFillingLevel(obj)           
            obj.f = (obj.m/obj.rho)/obj.V;               
            if obj.f > 0.9999
                obj.m = obj.m/obj.f;
                obj.f = 1;
            end
        end

    end
    

    methods(Static)

        function w = massfrac(conce)
            load variables/kinetic M
            m_ges = 0; %or obj.conce(5)*M_ROH
            for i =1:7
                m_ges = m_ges + conce(i)*M(i);
            end
            w = conce.*M./m_ges;        
        end

        function conce = setup_Reaction()
        %Setup for start concentration based on given ratio
            load variables/variables.mat  ratio
            conce = setup_cstart(ratio);
        end
   
    end
end