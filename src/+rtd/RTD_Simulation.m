%% Simulation model for reactive extrusion in twin screw extruder
% Sets up and simulates a twin screw extruder with user defined geometry
% and process conditions. The model uses a C-Chamber approach combined in
% an ideal CSTR cascade.
% Parameters for the ring opening polymerization of lactide to polylactide
% are stored.

tic
clear all
close all
clc
global e

%% User defined inputs for extruder setup and model
% set process conditions
%
% m_Flow0 = 50*1e-3/3600;    % input mass flow in kg/s 
m_Flow0 = 1e-4;
N = 100/60;         % screw rpm

% numerics related simulation inputs 
time = 10800;             % Simulation end time in seconds
plot_step = 5;          % plot only n-th time step (plotting every solution step is computationally very expensive)

% reaction related simulation inputs
T = 25;                 % Temperature in C
T_melt = 98;            % Melt temperature of lactide in C
P_U = 1e5;              % Ambient pressure in Pa
ratio = [500;2;1];      %    Molar ratio of monomer:co-initiator:catalyst
is_reaction = 'off';    % Set wheter reaction takes place ('on') or not ('off')

save('+rtd\+data\+variables\variables.mat', 'T', 'P_U', 'm_Flow0', 'N', 'ratio', 'is_reaction', 'T_melt');

% set extruder geometry
D_e = 11e-3;        % outer radius of extruder screw in m
D_i = 6e-3;         % inner radius of extruder screw in m
D_b = 11.5e-3;      % radius of barrel in m
C_l = 9e-3;         % distance between extruder screws in m

%% Screw configuration
% Example configuration for right handed screw element:
RSE1.type = 1;
RSE1.pitch = 1*D_e;         % pitch of screw element
RSE1.n_F = 2;               % number of flights
RSE1.n_reac = 20;            % number of reactors to model RSE
RSE1.f_ini = 0;             % initial filling level of RSE
RSE1.L = 10*D_e;             % length of RSE
RSE1.A1 = [];               % parameter to calculate drag flow, estimated when []
RSE1.A2 = 800;              % parameter to calculate pressure flow
RSE1.A3 = [];               % parameter to calculate shear rate, estimated when []

% Example configuration for kneading element KB90 (Not conveying element)
% KB2.type = 0;
% KB2.theta = 90;            % staggering angle
% KB2.n_D = 2;               % number of discs
% KB2.n_F = 2;               % number of flights
% KB2.n_reac = 4;            % number of reactors to model KB
% KB2.f_ini = 0;             % initial filling level of KB
% KB2.L = 1*D_e;             % length of KB
% KB2.A1 = 0;                % parameter to calculate drag flow
% KB2.A2 = 0;                % parameter to calculate pressure flow
% KB2.A0 = 5000;             % parameter to calculate pressure flow in neutral elements 
% KB2.A3 = [];               % parameter to calculate shear rate, estimated when []

% Example configuration for kneading element (Conveying element)
KB2.type = 3;
KB2.theta = 45;            % staggering angle
KB2.n_D = 2;               % number of discs
KB2.n_F = 2;               % number of flights
KB2.n_reac = 4;            % number of reactors to model KB
KB2.f_ini = 0;             % initial filling level of KB
KB2.L = 1*D_e;             % length of KB
KB2.A1 = 0.15;             % parameter to calculate drag flow
KB2.A2 = 250;              % parameter to calculate pressure flow

% Example configuration for left handed screw element:
LSE3.type = -1;
LSE3.pitch = 0.5*D_e;
LSE3.n_F = 2;               % number of flights
LSE3.n_reac = 4;            % number of reactors to model LSE
LSE3.f_ini = 0;             % initial filling level of LSE
LSE3.L = 1*D_e;             % length of LSE
LSE3.A1 = [];               % parameter to calculate drag flow, estimated when []
LSE3.A2 = 800;              % parameter to calculate pressure flow
LSE3.A3 = [];               % parameter to calculate shear rate, estimated when []

% Example configuration for right handed screw element:
RSE4.type = 1;
RSE4.pitch = 1*D_e;
RSE4.n_F = 2;               % number of flights
RSE4.n_reac = 10;            % number of reactors to model RSE
RSE4.f_ini = 0;             % initial filling level of RSE
RSE4.L = 5*D_e;             % length of RSE
RSE4.A1 = [];               % parameter to calculate drag flow, estimated when []
RSE4.A2 = 800;              % parameter to calculate pressure flow
RSE4.A3 = [];               % parameter to calculate shear rate, estimated when []

% % Example configuration for kneading element KB90 (Not conveying element)
% KB5.type = 0;
% KB5.theta = 90;            % staggering angle
% KB5.n_D = 2;               % number of discs
% KB5.n_F = 2;               % number of flights
% KB5.n_reac = 4;            % number of reactors to model KB
% KB5.f_ini = 0;             % initial filling level of KB
% KB5.L = 1*D_e;             % length of KB
% KB5.A1 = 0;                % parameter to calculate drag flow
% KB5.A2 = 0;                % parameter to calculate pressure flow
% KB5.A0 = 5000;             % parameter to calculate pressure flow in neutral elements 
% KB5.A3 = [];               % parameter to calculate shear rate, estimated when []

% Example configuration for kneading element KB45 (Conveying element)
KB5.type = 3;
KB5.theta = 45;            % staggering angle
KB5.n_D = 2;               % number of discs
KB5.n_F = 2;               % number of flights
KB5.n_reac = 4;            % number of reactors to model KB
KB5.f_ini = 0;             % initial filling level of KB
KB5.L = 1*D_e;             % length of KB
KB5.A1 = 0.15;             % parameter to calculate drag flow
KB5.A2 = 250;              % parameter to calculate pressure flow

% Example configuration for right handed screw element:
RSE6.type = 1;
RSE6.pitch = 1*D_e;
RSE6.n_F = 2;               % number of flights
RSE6.n_reac = 10;            % number of reactors to model RSE
RSE6.f_ini = 0;             % initial filling level of RSE
RSE6.L = 5*D_e;             % length of RSE
RSE6.A1 = [];               % parameter to calculate drag flow, estimated when []
RSE6.A2 = 800;              % parameter to calculate pressure flow
RSE6.A3 = [];               % parameter to calculate shear rate, estimated when []

% Example configuration for die (only cylindrical, one channel):
Die7.type = 2;
Die7.d_in = 5e-3;          % diameter of die at input
Die7.d_out = 5e-4;          % diameter of die at outout
Die7.n_reac = 4;            % number of reactors to model die
Die7.f_ini = 0;             % initial filling level of die
Die7.L = 2*D_e;             % length of die

% Arrange defined screw elements on screw
screw.E1 = RSE1;
screw.E2 = KB2;
screw.E3 = LSE3;
screw.E4 = RSE4;
screw.E5 = KB5;
screw.E6 = RSE6;
screw.E7 = Die7;

% Barrel temperature distribution
screw_temp = struct2cell(screw);
N_tot = 0;                                                  % total number of reactors
for i = 1:length(screw_temp), N_tot = N_tot + screw_temp{i}.n_reac; end
N_dist = N_tot - screw_temp{end}.n_reac;                    % number of reactor with equally distributed temperature along extruder
n = 7;                                                      % number of equally distributed temperature sections along extruder
T_barrel = [150 150 140 140 150 140 150 150];               % temperature distribution along extruder
T_b = rtd.helper.distributeVector(N_dist, n, T_barrel(1:end-1));
T_b(N_dist+1:N_tot) = T_barrel(end);


%% Set up extruder and simulation framework

% Define and create extruder 
e = rtd.Extruder(D_e,D_i,D_b,C_l, T_b);
e = designExtruder(e,screw,N);

% Simulate flow of extruder
m_0 = 0;
T_0 = T;
[t, m_sol] = rtd.CSTR_flow(e, m_0, T_0, [0 time], N_tot, plot_step);

% Simulate species
% load last_simulation.mat
% load kin.mat
c_in = [7.7277 0.015 0.003 0 0 0 0]; % Eingangskonzentration
[c, M_w, w_p] = rtd.CSTR_Reaction(e,c_in, [0 1000]);

% save last_simulation.mat
% load last_simulation.mat

% Simulate the extruder with temperature dependent viscosity
Reaction_data.c = c(end,:);
Reaction_data.M_w = M_w;
Reaction_data.w_p = w_p;
m_0 = m_sol(end,1:8:end);
T_0 = m_sol(end,8:8:end);
[t_2, m_sol_2] = rtd.CSTR_flow(e, m_0, T_0, [0 time], N_tot, plot_step, Reaction_data);


% Calculate mean residence time
filled_volume = 0;
for i = 1 : length(e.reactorConfiguration)
    filled_volume = filled_volume+e.reactorConfiguration{i}.V*e.reactorConfiguration{i}.f;
end
tau_mean = filled_volume/(e.m_Flow0/e.reactorConfiguration{1}.rho);
elapsed_time = toc;



% Saving inputs
projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
results_dir = fullfile(projectRoot, 'results');
if ~exist(results_dir, 'dir'); mkdir(results_dir); end
runID = char(datetime('now','Format','yyyyMMdd_HHmmss'));   % z.B. 20250813_101544
results_dir_actual = fullfile(results_dir, runID);
if ~exist(results_dir_actual, 'dir'); mkdir(results_dir_actual); end
save(fullfile(results_dir_actual, 'last_simulation.mat'), 'runID');
