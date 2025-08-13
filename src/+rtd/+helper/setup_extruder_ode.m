function dydt = setup_extruder_ode(~,y)
% setup the system of ordinary differtential equations according to the
% designed extruder e for the simulation of time dependent behavior
global e
% set starting values for e according to current solution of ode
for ixReactor = 1:length(e.reactorConfiguration)
    e.reactorConfiguration{ixReactor}.m = y((ixReactor-1)*8+1:(ixReactor-1)*8+7);
    e.reactorConfiguration{ixReactor}.T_m(2) = y((ixReactor-1)*8+8);
    e.reactorConfiguration{ixReactor} = calcFillingLevel(e.reactorConfiguration{ixReactor});
    e.reactorConfiguration{ixReactor}.w(:,2) = e.reactorConfiguration{ixReactor}.m/sum(e.reactorConfiguration{ixReactor}.m);
    e.reactorConfiguration{ixReactor}.w(isnan(e.reactorConfiguration{ixReactor}.w)) = 0;
    e.reactorConfiguration{ixReactor}.m_Flow = [0 0 0 0; 0 0 0 0];
    e.reactorConfiguration{ixReactor} = calcMolarMasses(e.reactorConfiguration{ixReactor});
end
e.reactorConfiguration{1}.m_Flow = [e.m_Flow0 0 0 0; 0 0 0 0];

% calculate drag flows for all reactors and add parameters to vector
k_p = zeros(length(e.reactorConfiguration),1);
eta_av = zeros(length(e.reactorConfiguration),1);
rho = zeros(length(e.reactorConfiguration),1);
L = zeros(length(e.reactorConfiguration),1);
D = zeros(length(e.reactorConfiguration),1);
for ixReactor = 1:length(e.reactorConfiguration)
    e.reactorConfiguration{ixReactor} = calcDFlow(e.reactorConfiguration{ixReactor});
    k_p(ixReactor) = e.reactorConfiguration{ixReactor}.k_p;
    e.reactorConfiguration{ixReactor} = calcViscosityPLA(e.reactorConfiguration{ixReactor});
    eta_av(ixReactor) = e.reactorConfiguration{ixReactor}.eta_av;    
    rho(ixReactor) = e.reactorConfiguration{ixReactor}.rho;
    L(ixReactor) = e.reactorConfiguration{ixReactor}.L;
    D(ixReactor) = e.reactorConfiguration{ixReactor}.D_e;
end

% connect streams between reactors
e.reactorConfiguration{1}.m_Flow(:,4) = e.reactorConfiguration{2}.m_Flow(:,2);
for i = 2:length(e.reactorConfiguration)
    e.reactorConfiguration{i}.m_Flow(:,1) = e.reactorConfiguration{i-1}.m_Flow(:,3);
    if i < length(e.reactorConfiguration)
        e.reactorConfiguration{i}.m_Flow(:,4) = e.reactorConfiguration{i+1}.m_Flow(:,2);
    end
end

A = zeros(length(e.reactorConfiguration));
B = zeros(length(e.reactorConfiguration),1);
A(1,1) = 1;
B(1) = e.P_U;
for i = 2:length(e.reactorConfiguration)
    if  e.reactorConfiguration{i}.f < 1
        A(i,i) = 1;
        B(i) = e.P_U;
    else
        A(i,i-1) = -rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i));
        A(i,i) = 2*rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i));
        A(i,i+1) = -rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i));
        B(i) = e.reactorConfiguration{i}.m_Flow(1,1)+e.reactorConfiguration{i}.m_Flow(1,4)...
            -e.reactorConfiguration{i}.m_Flow(1,3)-e.reactorConfiguration{i}.m_Flow(1,2);
    end
end
A(i+1,i+1) = 1;
B(i+1) = e.P_U;
p_sol = linsolve(A,B);
for i = 1:length(e.reactorConfiguration)
   e.reactorConfiguration{i}.p(2) = p_sol(i); 
   e.reactorConfiguration{i}.p(e.reactorConfiguration{i}.p < 1e5) = 1e5;
end

% connect pressure between reactors
e.reactorConfiguration{1}.p(3) = e.reactorConfiguration{2}.p(2);
e.reactorConfiguration{1}.p(1) = e.P_U;
for i = 2:length(e.reactorConfiguration)
    e.reactorConfiguration{i}.p(1) = e.reactorConfiguration{i-1}.p(2);
    if i < length(e.reactorConfiguration)
        e.reactorConfiguration{i}.p(3) = e.reactorConfiguration{i+1}.p(2);
    end
end

% calculate pressure profile with reactor specific parameters
A_new = zeros(length(e.reactorConfiguration));
A_new(1,1) = 1;
for i = 2:length(e.reactorConfiguration)
    if  e.reactorConfiguration{i}.f < 1
        A_new(i,i) = 1;
    else
        if e.reactorConfiguration{i}.p(2) > e.reactorConfiguration{i}.p(1)...
           && e.reactorConfiguration{i}.p(2) >= e.reactorConfiguration{i}.p(3)
                    A_new(i,i-1) = -rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i));
                    A_new(i,i) = 2*rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i));
                    A_new(i,i+1) = -rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i));
        elseif e.reactorConfiguration{i}.p(2) >= e.reactorConfiguration{i}.p(1)...
                && e.reactorConfiguration{i}.p(2) < e.reactorConfiguration{i}.p(3)
                    A_new(i,i-1) = -rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i));
                    A_new(i,i) = rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i))+rho(i+1)*D(i+1)^4/(k_p(i+1)*L(i+1)*eta_av(i+1));
                    A_new(i,i+1) = -rho(i+1)*D(i+1)^4/(k_p(i+1)*L(i+1)*eta_av(i+1));
        elseif e.reactorConfiguration{i}.p(2) < e.reactorConfiguration{i}.p(1)...
                && e.reactorConfiguration{i}.p(2) >= e.reactorConfiguration{i}.p(3)
                    A_new(i,i-1) = -rho(i-1)*D(i-1)^4/(k_p(i-1)*L(i-1)*eta_av(i-1));
                    A_new(i,i) = rho(i-1)*D(i-1)^4/(k_p(i-1)*L(i-1)*eta_av(i-1))+rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i));
                    A_new(i,i+1) = -rho(i)*D(i)^4/(k_p(i)*L(i)*eta_av(i));
        elseif e.reactorConfiguration{i}.p(2) < e.reactorConfiguration{i}.p(1)...
                && e.reactorConfiguration{i}.p(2) < e.reactorConfiguration{i}.p(3)
                    A_new(i,i-1) = -rho(i-1)*D(i-1)^4/(k_p(i-1)*L(i-1)*eta_av(i-1));
                    A_new(i,i) = rho(i-1)*D(i-1)^4/(k_p(i-1)*L(i-1)*eta_av(i-1))+rho(i+1)*D(i+1)^4/(k_p(i+1)*L(i+1)*eta_av(i+1));
                    A_new(i,i+1) = -rho(i+1)*D(i+1)^4/(k_p(i+1)*L(i+1)*eta_av(i+1));
        end
    end
end
A_new(i+1,i+1) = 1;
[p_sol,~] = linsolve(A_new,B);
for i = 1:length(e.reactorConfiguration)
   e.reactorConfiguration{i}.p(2) = p_sol(i); 
   e.reactorConfiguration{i}.p(e.reactorConfiguration{i}.p < 1e5) = 1e5;
end

% connect pressure between reactors
e.reactorConfiguration{1}.p(3) = e.reactorConfiguration{2}.p(2);
e.reactorConfiguration{1}.p(1) = e.P_U;
for i = 2:length(e.reactorConfiguration)
    e.reactorConfiguration{i}.p(1) = e.reactorConfiguration{i-1}.p(2);
    if i < length(e.reactorConfiguration)
        e.reactorConfiguration{i}.p(3) = e.reactorConfiguration{i+1}.p(2);
    end
end

% calculate pressure flows based on pressure profile
for i = 2:length(e.reactorConfiguration)
    if e.reactorConfiguration{i}.f == 1
        e.reactorConfiguration{i} = e.reactorConfiguration{i}.calcPFlow;
    end
end

% connect streams between reactors
e.reactorConfiguration{1}.m_Flow(:,4) = e.reactorConfiguration{2}.m_Flow(:,2);
e.reactorConfiguration{1}.w(:,3) = e.reactorConfiguration{2}.w(:,2);
for i = 2:length(e.reactorConfiguration)
    e.reactorConfiguration{i}.m_Flow(:,1) = e.reactorConfiguration{i-1}.m_Flow(:,3);
    e.reactorConfiguration{i}.w(:,1) = e.reactorConfiguration{i-1}.w(:,2);
    if i < length(e.reactorConfiguration)
        e.reactorConfiguration{i}.m_Flow(:,4) = e.reactorConfiguration{i+1}.m_Flow(:,2);
        e.reactorConfiguration{i}.w(:,3) = e.reactorConfiguration{i+1}.w(:,2);
    end
    % Adjust streams between reactors if mass flow balance is not kept.
    if (sum(e.reactorConfiguration{i}.m_Flow(:,1)) - sum(e.reactorConfiguration{i}.m_Flow(:,2)))>e.m_Flow0
         e.reactorConfiguration{i}.m_Flow(2,1) = e.m_Flow0 - e.reactorConfiguration{i}.m_Flow(1,1)...
             +sum(e.reactorConfiguration{i}.m_Flow(:,2));
    end
    if (sum(e.reactorConfiguration{i}.m_Flow(:,3)) - sum(e.reactorConfiguration{i}.m_Flow(:,4)))>e.m_Flow0
         e.reactorConfiguration{i}.m_Flow(2,3) = e.m_Flow0 - e.reactorConfiguration{i}.m_Flow(1,3)...
             +sum(e.reactorConfiguration{i}.m_Flow(:,4));
    end
end

% connect temperatures (T_m) between reactors
e.reactorConfiguration{1}.T_m(3) = e.reactorConfiguration{2}.T_m(2);
e.reactorConfiguration{1}.T_m(1) = e.T;
for i = 2:length(e.reactorConfiguration)
    e.reactorConfiguration{i}.T_m(1) = e.reactorConfiguration{i-1}.T_m(2);
    if i < length(e.reactorConfiguration)
        e.reactorConfiguration{i}.T_m(3) = e.reactorConfiguration{i+1}.T_m(2);
    end
end

% calculate heat flows for all reactors
for ixReactor = 1:length(e.reactorConfiguration)
    e.reactorConfiguration{ixReactor} = calcHeatFlow_bm(e.reactorConfiguration{ixReactor});
    e.reactorConfiguration{ixReactor} = calcheatFlow_diss(e.reactorConfiguration{ixReactor});
end


%% set areas in extruder where reaction should be on or off
% load variables/variables.mat T_melt
% for ixReactor = 2:length(e.reactorConfiguration)
%     % Initiallization
%     if (e.reactorConfiguration{ixReactor-1}.T_m(2)<T_melt &&...
%             e.reactorConfiguration{ixReactor}.T_m(2)>=T_melt &&...
%             strcmp(e.reactorConfiguration{ixReactor}.is_reaction, 'off'))
%         conce_start = Reaktor.setup_Reaction();
%         e.reactorConfiguration{ixReactor}.w(:,1) = Reaktor.massfrac(conce_start);
%         indexReaction = ixReactor;
%     end
%    % setting reaction in reactors "on" or "off"
%    if exist('indexReaction', 'var') && ixReactor>=indexReaction
%        e.reactorConfiguration{ixReactor}.is_reaction = 'on';
%    end
% end

%% define dydt
e.reactorConfiguration{1} = ode(e.reactorConfiguration{1});
dydt = e.reactorConfiguration{1}.dydt;
for ixR = 2:length(e.reactorConfiguration)
    e.reactorConfiguration{ixR} = ode(e.reactorConfiguration{ixR}); 
    dydt = [dydt; e.reactorConfiguration{ixR}.dydt];  
end
end