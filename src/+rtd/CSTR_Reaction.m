function [c, M_w, w_p] = CSTR_Reaction(extruder, c_in, tspan)

    reactorStruct = rtd.helper.reactor_struct_for_reaction(extruder);
    
    global c0;
    % Anzahl der Reaktoren und Komponenten
    n_reactors = numel(reactorStruct);
    n_species = numel(reactorStruct(1).InitialConcentration);
    
    % Initialbedingungen
    c0 = [];
    for i = 1:n_reactors
        c0 = [c0; reactorStruct(i).InitialConcentration(:)];
    end
    
    % Berechnung der Geschwindigkeitsraten (k)
    load +rtd\+data\+variables\variables.mat ratio
    load +rtd\+data\+variables\kinetic.mat C1_kin
    load +rtd\+data\+variables\kin.mat
    M0 = 7.7277e3;
    ROH0 = M0/ratio(1)*ratio(2);
    k_in = k;
    k = zeros(8, n_reactors);
    for i = 1:n_reactors
        if (extruder.reactorConfiguration{i}.T_m(2) > extruder.reactorConfiguration{i}.T_melt)
            k(:,i) = rtd.helper.setup_kinparameters(extruder.reactorConfiguration{i}.T_m(2), ROH0 ,C1_kin);
        else
            k(:,i) = k_in;
        end
       
    end
    
    % DGL-System lsen
    % options = odeset('RelTol',1e-6,'AbsTol',1e-8);
    options = odeset('OutputFcn',@(t,y,flag) myOutputFcn(t,y,flag),'RelTol',1e-6,'AbsTol',1e-8);
    [t, c] = ode15s(@(t,c) reactorSystem(t, c, reactorStruct, n_species, k, c_in), tspan, c0, options);      
    
    % Berechnung der Masse nach c
    c_last = c(end,:)';
    m = zeros(n_species, n_reactors);
    M_w = zeros(n_reactors,1);
    w_p = zeros(n_reactors,1);
    for i=1:n_reactors
        volume = extruder.reactorConfiguration{i}.V*extruder.reactorConfiguration{i}.f*1000;
        extruder.reactorConfiguration{i}.mol = extruder.reactorConfiguration{i}.m./ extruder.reactorConfiguration{i}.molarMass;
        extruder.reactorConfiguration{i} = calcMolarMasses(extruder.reactorConfiguration{i});
        m(:,i) = volume * c_last((i-1)*7+1:(i-1)*7+7) .* extruder.reactorConfiguration{i}.molarMass;
        extruder.reactorConfiguration{i}.mol = m(:,i)./ extruder.reactorConfiguration{i}.molarMass;
        extruder.reactorConfiguration{i} = calcMolarMasses(extruder.reactorConfiguration{i});
        m(:,i) = volume * c_last((i-1)*7+1:(i-1)*7+7) .* extruder.reactorConfiguration{i}.molarMass;
        M_w(i) = extruder.reactorConfiguration{i}.Mn(2) * 1.2;
        w_p(i) = (m(6,i)+m(7,i))/sum(m(:,i));
    end
    fprintf('Massenvektor fr Viskosittsberechnung: \n \n');
    disp(array2table(m));
end

function dc = reactorSystem(t, c, reactorStruct, n_species, k, c_in)
    n_reactors = numel(reactorStruct);
    dc = zeros(n_reactors*n_species, 1);
    
    for i = 1:n_reactors
        reactor = reactorStruct(i);
        idx = (i-1)*n_species+1 : i*n_species;
        c_current = c(idx);
        V = reactor.Volume*1000;
        
        inflow_terms = zeros(n_species,1);
        
        % Externer Zulauf
        if i == 1
            inflow_terms = inflow_terms + reactor.Inflow * c_in(:);
        end
        
        % Interne Zuflsse
        for j = 1:numel(reactor.Inflows)
            src = reactor.Inflows(j).Reactor;
            Q = reactor.Inflows(j).FlowRate;
            src_idx = (src-1)*n_species+1 : src*n_species;
            inflow_terms = inflow_terms + Q.* c(src_idx);
        end
        
        % Ausstrombehandlung
        % if i == n_reactors
        %     outflow_total = reactor.Outflows.FlowRate;
        % else
            outflow_total = sum([reactor.Outflows.FlowRate]);
        % end
        
        % Massenbilanz
        reaction_rates = rtd.helper.dcdt_reaction(c_current, k(:,i));
        dc(idx) = (inflow_terms - outflow_total*c_current)/V + reaction_rates;
    end
    t;
end

function status = myOutputFcn(t,y,flag)
    global c0
    persistent conversion
    switch flag  
        case 'init'
            figure
            xlabel('reactor')
            ylabel('conversion')
            conversion = animatedline;
        case []
            for i = 1:length(y)/7
                conv(i) = (c0(1)-y((i-1)*7+1))/c0(1);
            end
            clearpoints(conversion)
            addpoints(conversion,1:i,conv)
            drawnow
        case 'done'
    end
    status = 0;
end
