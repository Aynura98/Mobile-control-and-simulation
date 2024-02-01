%% parameters.m
% Beacon positions
S=  [0 0;
    14 0;
    14 14;
    0 14;
    8 6];

W = diag([0.01 0.01 deg2rad(5)].^2);
% W=diag([0.1 0.1 deg2rad(5)].^2);
V = eye(size(S,1))*0.01^2;
%V=eye(size(S,1))*0.5^2;

P0 = diag([0.15 0.15 deg2rad(10)].^2);
% P0=diag([2 2 deg2rad(40)].^2);

% Ts=1;
load path_planning_variables.mat nsteps T;
Ts = 0.75*(T/nsteps);       % Sampling time