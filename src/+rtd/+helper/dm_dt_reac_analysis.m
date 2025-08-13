clc; 
close all; 
% clear all;
% 
% global e
% 
load last_simulation_reaction.mat

% Reaktor zur Analyse (Eingabe)
reactor = input('Bitte Reaktor-Nummer eingeben (1-56): ');


%% Variablendeklaration
n_reac = 1:1:length(e.reactorConfiguration);
x = 1:n_reac(end);
x_diff = 1.5:n_reac(end);
L_extruder = zeros(size(n_reac));
L_extruder(1) = e.reactorConfiguration{1}.L;
for i = 2:n_reac(end)
    L_extruder(i) = L_extruder(i-1)+e.reactorConfiguration{i}.L;
end  
f_actual = zeros(length(e.reactorConfiguration),1);
for idx = 1:length(e.reactorConfiguration)
    f_actual(idx) = e.reactorConfiguration{idx}.f;
end
t = extra_data.time;
f = zeros(length(t),1);
dc_dt = zeros(length(t),7);
V = zeros(length(t),1);
M = zeros(length(t),7);
c = zeros(length(t),7);
dV_dt = zeros(length(t),1);
dM_dt = zeros(length(t),7);
dm_dt_reac = zeros(length(t),7);
for i = 1 : length(t)
    dm_dt_reac(i,:) = extra_data.reactor_data(i).dmdt_analysis(reactor,:);
    f(i) = extra_data.reactor_data(i).f(reactor);
    dc_dt(i,:) = extra_data.reactor_data(i).dc_dt(reactor,:);
    V(i) = extra_data.reactor_data(i).V(reactor);
    M(i,:) = extra_data.reactor_data(i).M(reactor,:);
    c(i,:) = extra_data.reactor_data(i).c(reactor,:);
    dV_dt(i) = extra_data.reactor_data(i).dV_dt(reactor);
    dM_dt(i,:) = extra_data.reactor_data(i).dM_dt(reactor,:);
end

term1 = dc_dt.*V.*M;
term2 = c.*dV_dt.*M;
term3 = zeros(length(t),7);
term3(:,6:7) = c(:,6:7).*V.*dM_dt(:,6:7);

dm_dt_reac = term1 + term2 + term3;

m_last = zeros(length(e.reactorConfiguration),1);
T_last = zeros(length(e.reactorConfiguration),1);
for idx = 1:length(e.reactorConfiguration)
    m_last(idx) = sum(e.reactorConfiguration{idx}.m);
    T_last(idx) = e.reactorConfiguration{idx}.T_m(2);
end

%% Plotten
figure(1)
plot(f_actual, 'b', 'LineWidth',2)
hold on
xline(reactor, 'k', 'LineWidth',2)
xlabel('Reaktor')
ylim([0, 1.1])
ylabel('Fllstand')


component = input('Bitte Komponenten-Nummer eingeben (1-7): ');


figure(2)
subplot(3,1,1)
plot(t, dm_dt_reac(:,1), 'r-', 'LineWidth',1)
hold on
plot(t, dm_dt_reac(:,2), 'g-', 'LineWidth',1)
plot(t, dm_dt_reac(:,3), 'b-', 'LineWidth',1)
plot(t, dm_dt_reac(:,4), 'c-', 'LineWidth',1)
plot(t, dm_dt_reac(:,5), 'm-', 'LineWidth',1)
plot(t, dm_dt_reac(:,6), 'y-', 'LineWidth',1)
plot(t, dm_dt_reac(:,7), 'k-', 'LineWidth',1)
plot(t, sum(dm_dt_reac(:,1:7),2), 'Color', [8 88 162]/255, 'LineWidth',2)
% set(gca, 'YScale', 'log')
grid
xlabel('Zeit / s')
ylabel('d(\cdot)/dt')
title(sprintf('dm/dt_{Reaktion}. Reaktor %d', reactor))
legend({'i=1 (M)','i=2 (ROH)','i=3 (C)','i=4 (I)','i=5 (A)','i=6 (R)','i=7 (D)', 'Summe(i)'}, ...
       'Orientation', 'horizontal', 'NumColumns', 4);


subplot(3,1,2)
plot(t, term1(:,component), 'r-', 'LineWidth',1)
hold on
plot(t, term2(:,component), 'g-', 'LineWidth',1)
plot(t, term3(:,component), 'b-', 'LineWidth',1)
% set(gca, 'YScale', 'log')
grid
xlabel('Zeit / s')
ylabel('d(\cdot)/dt')
title(sprintf('dm_{%d}/dt_{Reaktion} = term1_%d + term2_%d + term3_%d', component, component, component, component))
legend({'term1','term2','term3'},'Orientation', 'horizontal');

subplot(3,1,3)
plot(t, dc_dt(:,component), 'r-', 'LineWidth',1)
hold on
plot(t, c(:,component), 'c-', 'LineWidth',1)
plot(t, dV_dt, 'm-', 'LineWidth',1)
plot(t, V, 'g-', 'LineWidth',1)
plot(t, dM_dt(:,component), 'y-', 'LineWidth',1)
plot(t, M(:,component), 'b-', 'LineWidth',1)
% set(gca, 'YScale', 'log')
grid
xlabel('Zeit / s')
ylabel('d(\cdot)/dt')
title(sprintf('Detail Komponente: %d', component))
legend({'dc/dt','c','dV/dt','V', 'dM/dt', 'M'},'Orientation', 'horizontal');

figure(3)
subplot(2,1,1)
plot(x, m_last, 'r-', 'LineWidth',1)
hold on
plot(x_diff, diff(m_last), 'b-', 'LineWidth',1)
yline(0, 'k', 'LineWidth',1)
grid
xlabel('Reaktor')
ylabel('sum(m)')
legend('sum(m)','d(sum(m))/dx')
subplot(2,1,2)
plot(x, T_last, 'r-', 'LineWidth',1)
hold on
plot(x_diff, diff(T_last), 'b-', 'LineWidth',1)
yline(0, 'k', 'LineWidth',1)
grid
xlabel('Reaktor')
ylabel('T')
legend('T','d(T)/dx')

figure(4)
factor = max(diff(T_last))/max(diff(m_last));
plot((diff(T_last) + factor*diff(m_last)), 'b-', 'LineWidth',1)
yline(0, 'k', 'LineWidth',1)
grid
xlabel('Reaktor')
ylabel('d(y)/dx')
legend('d(y)/dx')