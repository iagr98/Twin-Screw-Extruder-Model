% set extruder geometry
D_e = 11e-3;        % outer radius of extruder screw in m
D_i = 6e-3;         % inner radius of extruder screw in m
D_b = 11.5e-3;      % radius of barrel in m
C_l = 9e-3;         % distance between extruder screws in m

%% Screw configuration
% 1. RSE (11 elements, 10mm pitch)
RSE1_pitch = 10e-3;
RSE1_n_elements = 11;
RSE1.type = 1;
RSE1.pitch = RSE1_pitch;                % pitch of screw element / m
RSE1.n_F = 2;                           % number of flights
RSE1.n_reac = 2*RSE1_n_elements;        % number of reactors to model RSE
RSE1.f_ini = 0;                         % initial filling level of RSE
RSE1.L = RSE1_n_elements*RSE1_pitch;    % length of RSE
RSE1.T_barrel = 160;                    % Barrel temperature
RSE1.A1 = [];                           % parameter to calculate drag flow, estimated when []
RSE1.A2 = 800;                          % parameter to calculate pressure flow
RSE1.A3 = [];                           % parameter to calculate shear rate, estimated when []

% 2. KB45 (10 elements,3mm length/disk)
KB2_disk_length = 3e-3;
KB2_n_elements = 10;
KB2.type = 3;
KB2.theta = 45;                         % staggering angle
KB2.n_D = KB2_n_elements;               % number of discs
KB2.n_F = 2;                            % number of flights
KB2.n_reac = 2*KB2.n_D;                 % number of reactors to model KB
KB2.f_ini = 0;                          % initial filling level of KB
KB2.L = KB2_n_elements*KB2_disk_length; % length of KB
KB2.T_barrel = 160;                     % Barrel temperature
KB2.A1 = 0.15;                          % parameter to calculate drag flow
KB2.A2 = 250;                           % parameter to calculate pressure flow

% 3. RSE (6 elements, 10mm pitch)
RSE3_pitch = 10e-3;
RSE3_n_elements = 6;
RSE3.type = 1;
RSE3.pitch = RSE3_pitch;                % pitch of screw element / m
RSE3.n_F = 2;                           % number of flights
RSE3.n_reac = 2*RSE3_n_elements;        % number of reactors to model RSE
RSE3.f_ini = 0;                         % initial filling level of RSE
RSE3.L = RSE3_n_elements*RSE3_pitch;    % length of RSE
RSE3.T_barrel = 160;                    % Barrel temperature
RSE3.A1 = [];                           % parameter to calculate drag flow, estimated when []
RSE3.A2 = 800;                          % parameter to calculate pressure flow
RSE3.A3 = [];                           % parameter to calculate shear rate, estimated when []

% 4. LSE (1 element, 5mm pitch)
LSE4_pitch = 10e-3/2;
LSE4_n_elements = 1;
LSE4.type = -1;
LSE4.pitch = LSE4_pitch;                % pitch of screw element / m
LSE4.n_F = 2;                           % number of flights
LSE4.n_reac = 4*LSE4_n_elements;        % number of reactors to model LSE
LSE4.f_ini = 0;                         % initial filling level of LSE
LSE4.L = 2*LSE4_n_elements*LSE4_pitch;  % length of LSE
LSE4.T_barrel = 160;                    % Barrel temperature
LSE4.A1 = [];                           % parameter to calculate drag flow, estimated when []
LSE4.A2 = 800;                          % parameter to calculate pressure flow
LSE4.A3 = [];                           % parameter to calculate shear rate, estimated when []

% 5. KB45 (6 elements, 3mm length/disk)
KB5_disk_length = 3e-3;
KB5_n_elements = 6;
KB5.type = 3;
KB5.theta = 45;                         % staggering angle
KB5.n_D = KB5_n_elements;               % number of discs
KB5.n_F = 2;                            % number of flights
KB5.n_reac = 2*KB5.n_D;                 % number of reactors to model KB
KB5.f_ini = 0;                          % initial filling level of KB
KB5.L = KB5_n_elements*KB5_disk_length; % length of KB
KB5.T_barrel = 160;                     % Barrel temperature
KB5.A1 = 0.15;                          % parameter to calculate drag flow
KB5.A2 = 250;                           % parameter to calculate pressure flow

% 6. RSE (6 elements, 10mm pitch)
RSE6_pitch = 10e-3;
RSE6_n_elements = 6;
RSE6.type = 1;
RSE6.pitch = RSE6_pitch;                % pitch of screw element / m
RSE6.n_F = 2;                           % number of flights
RSE6.n_reac = 2*RSE6_n_elements;        % number of reactors to model RSE
RSE6.f_ini = 0;                         % initial filling level of RSE
RSE6.L = RSE6_n_elements*RSE6_pitch;    % length of RSE
RSE6.T_barrel = 160;                    % Barrel temperature
RSE6.A1 = [];                           % parameter to calculate drag flow, estimated when []
RSE6.A2 = 800;                          % parameter to calculate pressure flow
RSE6.A3 = [];                           % parameter to calculate shear rate, estimated when []

% 7. LSE (1 element, 5mm pitch)
LSE7_pitch = 10e-3/2;
LSE7_n_elements = 1;
LSE7.type = -1;
LSE7.pitch = LSE7_pitch;                % pitch of screw element / m
LSE7.n_F = 2;                           % number of flights
LSE7.n_reac = 4*LSE7_n_elements;        % number of reactors to model LSE
LSE7.f_ini = 0;                         % initial filling level of LSE
LSE7.L = 2*LSE7_n_elements*LSE7_pitch;  % length of LSE
LSE7.T_barrel = 160;                    % Barrel temperature
LSE7.A1 = [];                           % parameter to calculate drag flow, estimated when []
LSE7.A2 = 800;                          % parameter to calculate pressure flow
LSE7.A3 = [];                           % parameter to calculate shear rate, estimated when []

% 8. KB45 (3 elements, 3mm length/disk)
KB8_disk_length = 3e-3;
KB8_n_elements = 3;
KB8.type = 3;
KB8.theta = 45;                         % staggering angle
KB8.n_D = KB8_n_elements;               % number of discs
KB8.n_F = 2;                            % number of flights
KB8.n_reac = 2*KB8.n_D;                 % number of reactors to model KB
KB8.f_ini = 0;                          % initial filling level of KB
KB8.L = KB8_n_elements*KB8_disk_length; % length of KB
KB8.T_barrel = 160;                     % Barrel temperature
KB8.A1 = 0.15;                          % parameter to calculate drag flow
KB8.A2 = 250;                           % parameter to calculate pressure flow

% 9. KB90 (9 elements, 3mm length/disk)
KB9_disk_length = 3e-3;
KB9_n_elements = 9;
KB9.type = 0;
KB9.theta = 90;                         % staggering angle
KB9.n_D = KB9_n_elements;               % number of discs
KB9.n_F = 2;                            % number of flights
KB9.n_reac = 2*KB9.n_D;                 % number of reactors to model KB
KB9.f_ini = 0;                          % initial filling level of KB
KB9.L = KB9_n_elements*KB9_disk_length; % length of KB
KB9.T_barrel = 160;                     % Barrel temperature
KB9.A1 = 0;                             % parameter to calculate drag flow
KB9.A2 = 0;                             % parameter to calculate pressure flow
KB9.A0 = 5000;                          % parameter to calculate pressure flow in neutral elements 
KB9.A3 = [];                            % parameter to calculate shear rate, estimated when []

% 10. RSE (9 elements, 10mm pitch)
RSE10_pitch = 10e-3;
RSE10_n_elements = 10;
RSE10.type = 1;
RSE10.pitch = RSE10_pitch;               % pitch of screw element / m
RSE10.n_F = 2;                           % number of flights
RSE10.n_reac = 2*RSE10_n_elements;       % number of reactors to model RSE
RSE10.f_ini = 0;                         % initial filling level of RSE
RSE10.L = RSE10_n_elements*RSE10_pitch;  % length of RSE
RSE10.T_barrel = 160;                    % Barrel temperature
RSE10.A1 = [];                           % parameter to calculate drag flow, estimated when []
RSE10.A2 = 800;                          % parameter to calculate pressure flow
RSE10.A3 = [];                           % parameter to calculate shear rate, estimated when []

% 11 Die
Die11.type = 2;
Die11.d_in = 5e-3;                      % diameter of die at input
Die11.d_out = 5e-4;                     % diameter of die at outout
Die11.n_reac = 4;                       % number of reactors to model die
Die11.f_ini = 0;                        % initial filling level of die
Die11.L = 2*D_e;                        % length of die
Die11.T_barrel = 160;                   % Barrel temperature


% Arrange defined screw elements on screw
screw.E1 = RSE1;
screw.E2 = KB2;
screw.E3 = RSE3;
screw.E4 = LSE4;
screw.E5 = KB5;
screw.E6 = RSE6;
screw.E7 = LSE7;
screw.E8 = KB8;
screw.E9 = KB9;
screw.E10 = RSE10;
screw.E11 = Die11;


save('+rtd\+data\screwGeometry.mat', 'screw', 'D_e','D_i', 'D_b', 'C_l')