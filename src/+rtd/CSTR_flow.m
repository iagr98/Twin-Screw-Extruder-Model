function [t, m] = CSTR_flow(extruder, m_ini, T_ini, t_span, N_tot, plot_step, Reaction_data)

global e
e = extruder;
T_melt = 98;


% Setup reaction start conditions
if nargin < 7 % checks if concentrations c were given as input
    plot_extruder = true;
else    
    c = Reaction_data.c;
    M_w = Reaction_data.M_w;
    w_p = Reaction_data.w_p;
    for i = 1 : length(e.reactorConfiguration)    % Berechnung und einsetzung der Viskositt
        if (e.reactorConfiguration{i}.T_m(2)>T_melt)
            e.reactorConfiguration{i}.M_w = M_w(i);                
            e.reactorConfiguration{i}.w_p = w_p(i);              
            e.reactorConfiguration{i} = calcViscosityPLA(e.reactorConfiguration{i});            
            e.reactorConfiguration{i}.c(:,2) = c((i-1)*7+1 : 7*i);
        end
    end
    for i = 1 : length(e.reactorConfiguration)    % Berechnung und einsetzung des reaktiven Terms (term4)
        if (e.reactorConfiguration{i}.T_m(2)>T_melt)            
            if i >= 2
                e.reactorConfiguration{i}.c(:,1) = e.reactorConfiguration{i-1}.c(:,2);
                if i < length(e.reactorConfiguration)
                    e.reactorConfiguration{i}.c(:,3) = e.reactorConfiguration{i+1}.c(:,2);
                end
            else
                e.reactorConfiguration{1}.c(:,3) = e.reactorConfiguration{2}.c(:,2);
            end     
            e.reactorConfiguration{i} = calcTerm4(e.reactorConfiguration{i}); 
            e.reactorConfiguration{i}.flag_term4 = 1;
        end        
        e.reactorConfiguration{i}.m = m_ini(i); 
    plot_extruder = false;
    end
end
 
% Set initial condition vector (y0) according to f_ini in screw and T
y0 = zeros(N_tot*2,1);
y0(1:2:end) = m_ini;
y0(2:2:end) = T_ini;

% Simulate the set up extruder
options = odeset('OutputFcn',@(t,y,flag) OutFcn(t,y,flag,plot_step, plot_extruder),...
    'NonNegative',1:length(y0),'RelTol',1e-5,'AbsTol',1e-7,'Events',@stationaryEvent);

[t,m_sol] = ode113(@(t,y) rtd.helper.setup_extruder_ode(t,y), t_span, y0, options);
m = m_sol;

end

%% outputfcn to plot extruder state during the simulation
function status = OutFcn(t,~,flag,plot_step, plot_extruder)
global e

persistent plot_reacts fill conversion pressure filling temperature plot_terms...
    terms L_extruder step_number viscosity heat_flow
  
switch flag
    case 'init'
        step_number = 1;
        % set up plots for filling during simulation
        figure
        plot_reacts = [1;2;3;4;5;6;7];
        colors = ["r";"g";"b";"c";"m";"y";"k"];
        legend_labels = "reactor" + " " + num2str(plot_reacts(1));
        subplot(3,2,2)
        for ip = 1:length(plot_reacts)
            legend_labels(ip) = "reactor"+" "+num2str(plot_reacts(ip));
            fill{ip} = animatedline('color',colors(ip));
        end
        xlabel('time / s')
        ylabel('fill level')
        legend(legend_labels)
        subplot(3,2,1)
        heat_flow = animatedline('color','k');        
        ylabel('Heat flow / J/s')
        xlim([0 inf])
        grid
        subplot(3,2,3)
        filling = animatedline('color','b');
        xlabel('length / m')
        ylabel('fill level / -')
        grid
        xlim([0 inf])
        subplot(3,2,4)
        pressure = animatedline;
        xlabel('length / m')
        ylabel('pressure / Pa')
        grid
        xlim([0 inf])
        subplot(3,2,5)
        temperature{1} = animatedline('color','r');
        temperature{2} = animatedline('color','k');
        xlabel('length / m')
        ylabel('T / C')
        grid
        xlim([0 inf])        
        subplot(3,2,6)
        viscosity = animatedline('color',[1,0.5,0]);
        xlabel('length / m')
        ylabel('\eta_{i,av}')
        grid
        xlim([0 inf])        
        % plot_terms = ["m(M)";"m(ROH)";"m(C)";"m(I)";"m(A)";"m(R)";"m(D)";"T"];
        % colors_terms = ["k";"b";"r";"g";"m";"y";"c";"b"];
        % legend_labels_terms = "Term" + " " + num2str(plot_terms(1));
        % for ip = 1:(length(plot_terms))
        %     legend_labels_terms(ip) = "\cdot ="+" "+num2str(plot_terms(ip));
        %     terms{ip} = animatedline('color',colors_terms(ip));
        % end
        % xlabel('length / m')
        % ylabel('m_{vor} m_{nach}')
        % ylabel('d(\cdot)/dt')
        % legend(legend_labels_terms)
        % grid   
        % xlim([0 inf])

        % Plot user defined extruder configuration
        n_reac = 1:1:length(e.reactorConfiguration);
        L_extruder = zeros(size(n_reac));
        L_extruder(1) = e.reactorConfiguration{1}.L;
        for i = 2:n_reac(end)
            L_extruder(i) = L_extruder(i-1)+e.reactorConfiguration{i}.L;
        end  
        if (plot_extruder)
            figure
            set(gca,'DataAspectRatio',[1 .01 1],'YTick',[])
            xlabel('dimensionless extruder length / L/D')
            title('screw configuration of simulated extruder');
            L_elements = zeros(length(e.screw),1);
            pitch_elements = zeros(length(e.screw),1);
            n_elements = zeros(length(e.screw),1);
            screw_type = zeros(length(e.screw),1);
            for i = 1:length(L_elements)
                L_elements(i) = e.screw{i,1}.L;
                screw_type(i) = e.screw{i,1}.type;
                switch screw_type(i)
                    case 1
                        pitch_elements(i) = e.screw{i,1}.pitch;
                        n_elements(i) = L_elements(i)/pitch_elements(i); 
                        for j = 1:n_elements(i)
                            x_end = sum(L_elements(1:i))-(n_elements(i)-j)*pitch_elements(i);
                            x_start = sum(L_elements(1:i))-(n_elements(i)-j+1)*pitch_elements(i);
                            x = [x_start,x_start,x_end,x_start];
                            y = [0,.05,0,0];
                            hold on
                            plot(x./e.extruder_geometry.D_e,y,'Color','b')
                        end
                    case -1
                        pitch_elements(i) = e.screw{i,1}.pitch;
                        n_elements(i) = L_elements(i)/pitch_elements(i); 
                        for j = 1:n_elements(i)
                            x_start = sum(L_elements(1:i))-(n_elements(i)-j)*pitch_elements(i);
                            x_end = sum(L_elements(1:i))-(n_elements(i)-j+1)*pitch_elements(i);
                            x = [x_start,x_start,x_end,x_start];
                            y = [0,.05,0,0];
                            hold on
                            plot(x./e.extruder_geometry.D_e,y,'Color','g')
                        end
                    case {0,3}
                        n_elements(i) = e.screw{i,1}.n_D; 
                        for j = 1:n_elements(i)
                            x_start = sum(L_elements(1:i-1))+(j-1)*L_elements(i)/n_elements(i);
                            x_end = sum(L_elements(1:i-1))+(j)*L_elements(i)/n_elements(i);
                            x = [x_start,x_start,x_end,x_end];
                            y = [0,.05,0.05,0];
                            hold on
                            plot(x./e.extruder_geometry.D_e,y,'Color','r')
                        end
                    case 2
                        x_start = sum(L_elements(1:i-1));
                        x_end = sum(L_elements(1:i));
                        x = [x_start,x_start,x_end,x_end];
                        y = [0,.05,0.01,0];
                        hold on
                        plot(x./e.extruder_geometry.D_e,y,'Color','black')
                end
            end
            xlim([0,L_extruder(end)/e.extruder_geometry.D_e])
        end
    
        % text(0,.06,num2str(1))
        % title('Extruder Configuration')
               
    case []
        if mod(step_number,plot_step) ~= 0 % draw solution only for every n-th step
           % don't draw anything 
        else
            vec_reac = 1:1:length(e.reactorConfiguration);
            for i = 1:length(e.reactorConfiguration)
                % store current extruder state in vectors
                p(i) = e.reactorConfiguration{i}.p(2);
                f(i) = e.reactorConfiguration{i}.f;
                conv(i) = e.reactorConfiguration{i}.conv;
                T_m(i) = e.reactorConfiguration{i}.T_m(2);             
                T_b(i) = e.reactorConfiguration{i}.T_b;             
                
                 % term0(i) = e.reactorConfiguration{i}.term0;
                 % term1(i) = e.reactorConfiguration{i}.term1;
                 % term2(i) = e.reactorConfiguration{i}.term2;
                 % term3(i) = e.reactorConfiguration{i}.term3;
                 % total(i) = term0(i) + term1(i) + term2(i) - term3(i);
                 % 
                 % mf0(i) = sum(e.reactorConfiguration{i}.m_Flow(:,1));
                 % mf1(i) = sum(e.reactorConfiguration{i}.m_Flow(:,2));
                 % mf2(i) = sum(e.reactorConfiguration{i}.m_Flow(:,3));
                 % mf3(i) = sum(e.reactorConfiguration{i}.m_Flow(:,4));
                 % m(i) = e.reactorConfiguration{i}.m(1);

                % dydt1(i) = e.reactorConfiguration{i}.dydt(1);
                % dydt2(i) = e.reactorConfiguration{i}.dydt(2);
                % dydt3(i) = e.reactorConfiguration{i}.dydt(3);
                % dydt4(i) = e.reactorConfiguration{i}.dydt(4);
                % dydt5(i) = e.reactorConfiguration{i}.dydt(5);
                % dydt6(i) = e.reactorConfiguration{i}.dydt(6);
                % dydt7(i) = e.reactorConfiguration{i}.dydt(7);
                % dydt8(i) = e.reactorConfiguration{i}.dydt(8);

                % m_forward(i) = sum(e.reactorConfiguration{i}.m_Flow(:,1)) - sum(e.reactorConfiguration{i}.m_Flow(:,2));
                % m_forward_2(i) = sum(e.reactorConfiguration{i}.m_Flow(:,3)) - sum(e.reactorConfiguration{i}.m_Flow(:,4));

                eta_av(i) = e.reactorConfiguration{i}.eta_av;
                % gamma_av(i) = e.reactorConfiguration{i}.gamma_av;
                dTdt(i) = e.reactorConfiguration{i}.dTdt;
                
                                                
            end
            for ip = 1:length(plot_reacts)
                % add points to the created plots 
                addpoints(fill{ip},t(end),e.reactorConfiguration{plot_reacts(ip)}.f)
            end
            clearpoints(heat_flow)
            addpoints(heat_flow,L_extruder,dTdt)
            clearpoints(filling)
            addpoints(filling,L_extruder,f)
            clearpoints(pressure)
            addpoints(pressure,L_extruder,p)            
            clearpoints(temperature{1})
            addpoints(temperature{1},L_extruder,T_m)   
            clearpoints(temperature{2})
            addpoints(temperature{2},L_extruder,T_b)
            clearpoints(viscosity)
            addpoints(viscosity,L_extruder,eta_av)

            
            % clearpoints(terms{1})            
            % addpoints(terms{1},L_extruder,dydt1)                        
            % clearpoints(terms{2})            
            % addpoints(terms{2},L_extruder,dydt2)                        
            % clearpoints(terms{3})
            % addpoints(terms{3},L_extruder,dydt3)
            % clearpoints(terms{4})
            % addpoints(terms{4},L_extruder,dydt4)                        
            % clearpoints(terms{5})
            % addpoints(terms{5},L_extruder,dydt5)                        
            % clearpoints(terms{6})
            % addpoints(terms{6},L_extruder,dydt6)                        
            % clearpoints(terms{7})
            % addpoints(terms{7},L_extruder,dydt7)                        
            % clearpoints(terms{8})
            % addpoints(terms{8},L_extruder,dydt8)                        

            drawnow
            if e.reactorConfiguration{1}.f == 1
                error('Error: First screw element full. Screw configuration not suitable for input mass flow.')
            end
        end
        step_number = step_number+1;
    case 'done'
end
status = 0;
end

%% Event function to stop simulation when stationary 
function [value,isterminal,direction] = stationaryEvent(t,y)
    persistent y_prev t_prev step_number % save solution of previous step
    if t == 0
        y_prev = y;
        step_number = 1;
    end
    if mod(step_number,100) ~= 0 % only check for every 1000th step
        value = 1; % Geben Sie 1 zurck, um die Ereignisfunktion fortzusetzen
        isterminal = 0; % Die Integration soll nicht gestoppt werden
        direction = 0; % Alle Richtungen zhlen
    else
        if isempty(y_prev) == 1 | isempty(t_prev) == 1
            y_prev = y; % Setzen Sie den vorherigen Zustand zu Beginn der Integration
            t_prev = t;
            value = 1; % Geben Sie 1 zurck, um die Ereignisfunktion fortzusetzen
            isterminal = 0; % Die Integration soll nicht gestoppt werden
            direction = 0; % Alle Richtungen zhlen
        else
            % Berechnen Sie die nderung der Lsung ber die Zeit
            change = (y-y_prev)./(t-t_prev); 
            if isnan(change)
                value = 1; % Geben Sie 1 zurck, um die Ereignisfunktion fortzusetzen
                isterminal = 1; % Stoppen Sie die Integration, wenn das Ereignis eintritt
                direction = 0; % Alle Richtungen zhlen        
                y_prev = y; % Aktualisieren Sie den vorherigen Zustand von y fr den nchsten Schritt
                t_prev = t;
            else
            
                % Schwer einen geeigneten Wert zu finden
                tolerance = 1e-8; % Toleranzschwelle fr die Stationaritt
        
                if all(abs(change)< tolerance)
                    value = 0; % Wert, der auslst, wenn das Ereignis eintritt
                else
                    value = 1; % Wert, der auslst, wenn das Ereignis nicht eintritt
                end
        
                isterminal = 1; % Stoppen Sie die Integration, wenn das Ereignis eintritt
                direction = 0; % Alle Richtungen zhlen
        
                y_prev = y; % Aktualisieren Sie den vorherigen Zustand von y fr den nchsten Schritt
                t_prev = t;
            end
        end
    end
    step_number = step_number+1;
end             