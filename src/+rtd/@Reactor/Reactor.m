classdef Reactor

    properties
        % Prozessbezogene Gren
        
        f
        m = 0
        T 
        T_b                 % Temperature of barrel
        T_m = [0 0 0]       % Temperature of the melt [T_m(i-1) T_m(i) T_m(i+1)]
        T_melt              % Melting temperature of polymer [C]
        gamma = 0
        lambda = 0.2        % Heat conductivity
        rho
        eta   
        eta_new
        c_p                 % heat capacity of polymer [J/kgK]
        mol = zeros(7,1)    % [Monomer, Co-Initiator, Catalyst, Initiator, Acid, R, D] / mol
        M0                  % Initial monomer concentration / mol/m3
        M_rep = 72e-3;      % molar mass of the monomeric repeating unit / kg/mol
        time = 0
        molarMass = zeros(7,1)
        c = zeros(7,3)      % Concentration of substance / mol/L
        Mn = zeros(3,1)
        N                   % Drehzal
        k                   % kinetic parameters of reaction
        ratio
        eta_params          % parameter for calculation of the viscosity
        number
        n_discret
        L                   % length of reactor
        V                   % maximum volume of reactor
        htc_b               % Heat transfer coefficient between barrel-melt
        HeatFlow_bm         % Heat flow between barrel and melt
        heatFlow_diss       % Heat flow per unit volume due to viscous dissipation
        melt_length         % Length of melting area in m
        
        
        C1_kin
        ROH0
        M_w = 10
        w_p = 1
        eta_min = 0.1

                
    end
    methods (Static)
        %Konstruktor der Klasse
        function obj = Reactor(obj)
            obj.rho = 1120;
            obj.c_p = 2260;    
            load +rtd\+data\+variables\variables.mat T T_melt M0
            load +rtd\+data\+variables\kinetic.mat M C1_kin ratio
            load +rtd\+data\+variables\viscosity.mat
            obj.eta_params.B = B;
            obj.eta_params.E_A = E_a;
            obj.eta_params.R = R;
            obj.eta_params.a = a;
            obj.eta_params.eta_0_ref = eta_0_ref;
            obj.ratio = ratio;
            obj.M0 = M0;
            obj.ROH0 = obj.M0/obj.ratio(1)*obj.ratio(2);
            obj.molarMass = M;
            obj.T = T;
            obj.T_melt = T_melt; 
            obj.C1_kin = C1_kin;
            %obj.k = rtd.helper.setup_kinparameters(150,obj.ROH0,obj.C1_kin);
                    
        end

        function [Mn, molarMass] = calcMolarMasses(obj, m, mol, molarMass)
            % calculates the molar masses of all polymer species according
            % to the current filling of the extruder reactor
            % molar masses of R and D are set to number average molar mass of polymer, may be changed when
            % moment equations are intoduced.
            volume = sum(m)/obj.rho;            
            conv = (obj.M0*volume-mol(1))/(obj.M0*volume);
            conv(conv<0)=0;
            molarMass(6:7) = (conv*obj.M0*volume)/((obj.M0/obj.ratio(1)*obj.ratio(2))*volume)*2*obj.M_rep+molarMass(4);
            molarMass(isnan(molarMass)) = molarMass(4);
            molarMass(isinf(obj.molarMass)) = molarMass(4);
            Mn = molarMass(6); % in kg/mol
        end

        function conce = setup_Reaction()
        %Setup for start concentration based on given ratio
            load variables/variables.mat  ratio
            conce = setup_cstart(ratio);
        end

        function k = arrhenius(Ea,k0,T)
            R = 8.314;
            k = k0*exp(-Ea/(R*(T+273.15)));
        end

        function M_eq = calc_Meq(T)
            dHlc = -23.3e3;
            dSlc = -22;
        %     dHlc = -50e3;
        %     dSlc = -12;
            M_eq = exp(dHlc/(8.314*(T+273.15))-dSlc/8.314);
        end

        function k = setup_kinparameters(T,ROH,kin)
        % calculates the kinetic parameters of the reaction system for a given
        % temperature T in C
            k_MC = rtd.Reactor.arrhenius(kin.kCM(1),kin.kCM(2),T);
            k_MC = k_MC + (kin.a) * ROH;
            % k_MC = 0;
            K_A = rtd.Reactor.arrhenius(kin.KA(1),kin.KA(2),T);
            k_A1 = 10e3;  
            % k_p = arrhenius(kin.kp(1),kin.kp(2),T);
            k_p = 300e3;
            M_eq = rtd.Reactor.calc_Meq(T);
            k_A2 = k_A1/K_A;
            k_d = k_p*M_eq;
            k_s = 1e6;
            k_te = rtd.Reactor.arrhenius(kin.kte(1),kin.kte(2),T);
            k_de = rtd.Reactor.arrhenius(kin.kde(1),kin.kde(2),T);
            k = [k_MC; k_A1; k_A2; k_p; k_d; k_s; k_te; k_de];
        end

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

        function reactor_test = reactor_struct_for_reaction(e)
            reactor_test = [];
            reactor_test(1).Volume = e.reactorConfiguration{1}.V;
            reactor_test(1).Inflows = [];
            reactor_test(1).Inflow = e.reactorConfiguration{1}.m_Flow(1,1);
            reactor_test(1).Outflows = [struct('Reactor', e.reactorConfiguration{1}.number+1, 'FlowRate',e.reactorConfiguration{1}.m_Flow(1,3)+e.reactorConfiguration{1}.m_Flow(2,3))];
            reactor_test(1).T = e.reactorConfiguration{1}.T_m(2);
            reactor_test(1).f = e.reactorConfiguration{1}.f;
            reactor_test(1).m = e.reactorConfiguration{1}.m;
            reactor_test(1).molarMass = e.reactorConfiguration{1}.molarMass;
            for i = 2:length(e.reactorConfiguration)-1
                reactor_test(i).Volume = e.reactorConfiguration{i}.V;
                reactor_test(i).Inflows = [struct('Reactor', e.reactorConfiguration{i}.number-1, 'FlowRate',e.reactorConfiguration{i}.m_Flow(1,1)+e.reactorConfiguration{i}.m_Flow(2,1)),...
                    struct('Reactor', e.reactorConfiguration{i}.number+1, 'FlowRate',e.reactorConfiguration{i}.m_Flow(1,4)+e.reactorConfiguration{i}.m_Flow(2,4))];
                reactor_test(i).Outflows = [struct('Reactor', e.reactorConfiguration{i}.number-1, 'FlowRate',e.reactorConfiguration{i}.m_Flow(1,2)+e.reactorConfiguration{i}.m_Flow(2,2)),...
                    struct('Reactor', e.reactorConfiguration{i}.number+1, 'FlowRate',e.reactorConfiguration{i}.m_Flow(1,3)+e.reactorConfiguration{i}.m_Flow(2,3))];
                reactor_test(i).T = e.reactorConfiguration{i}.T_m(2);
                reactor_test(i).f = e.reactorConfiguration{i}.f;
                reactor_test(i).m = e.reactorConfiguration{i}.m;
                reactor_test(i).molarMass = e.reactorConfiguration{i}.molarMass;
            end
            reactor_test(i+1).Volume = e.reactorConfiguration{i}.V;
            reactor_test(i+1).Inflows = [struct('Reactor', e.reactorConfiguration{i}.number, 'FlowRate',e.reactorConfiguration{i+1}.m_Flow(1,1)+e.reactorConfiguration{i+1}.m_Flow(2,1))];
            reactor_test(i+1).Outflows = [struct('Reactor', e.reactorConfiguration{i}.number, 'FlowRate',e.reactorConfiguration{i+1}.m_Flow(1,2)+e.reactorConfiguration{i+1}.m_Flow(2,2)),...
                    struct('Reactor', 0, 'FlowRate',e.reactorConfiguration{i+1}.m_Flow(1,3)+e.reactorConfiguration{i+1}.m_Flow(2,3))];
            reactor_test(i+1).T = e.reactorConfiguration{i}.T_m(2);
            reactor_test(i+1).f = e.reactorConfiguration{i}.f;
            reactor_test(i+1).m = e.reactorConfiguration{i}.m;
            reactor_test(i+1).molarMass = e.reactorConfiguration{i}.molarMass;
        end
   
    end
end