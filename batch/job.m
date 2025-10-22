clc, clear, close all
%% Add path to data
addpath(genpath('/home/zachariand/Exjobb'));
addpath(genpath('/home/zachariand/Exjobb/functions'));
%% Run Job
% FEM parameters

params.le = 0.0005; %halverad
params.lx = 0.1;
params.ly = 0.1;
params.wx = 0.04;
params.wy = 0.04;
params.loadcase = 3;
params.Vf = 0.3;

params.t = 0.1;
params.ngp = 4;

params.stressCon = 1;
params.pnm = 8;
% params.sigc = 1.15; %360e6; % Max stress for constraint
params.sigc = 1e-3;
params.ngr = 1;

params.R1tol = 1e-2;
% params.N = 6;
params.N = 1;
params.disp = -9e-4; % displacement [nodes total-size]

% Material parameters
E = 210e9;
v = 0.3;
params.E1 = E; params.E2 = 0.7*E; params.E3 = 0.5*E;
params.v12 = v; params.v13 = 0.5*v; params.v23 = 0.5*v;
% params.v21 = v; params.v32 = v; params.v31 = v;

params.sigy01 = 360e6;
params.sigy02 = 300e6;
params.sigy03 = 250e6;

params.H = 10e9;
params.Kinf = 0; %0.3*params.sigy01; %extra terms linHard (sat. stress)
params.xi = 1e-3; %extra terms linHard (sat. exp)

params.rtol = 1e-4;
params.PT = 1; % 1 for Deformation plasticity, 0 for Elastoplasticity

% Optimization Parameters
params.re = 6; % Elements in radius (dubblerad)
params.filtOn = true;
params.pad = true;
params.p = 1.5; %3;
params.q = 1; %2,5;
params.eta = 0.5;
params.beta = 0.1;
params.rampB = [2, 10, 1.13]; % [0/1/2, end value, factor size]  (0: off, 1: on, 2: on after simp converges)
params.rampPQ = [4, 0.1+2/30]; % [end value of p, increment size]
params.del = 1e-9;
params.dels = 1e-3;
%params.ncon = 1; % Nr of constraints
params.xtol = 1e-3;
params.iterMax = 1200;
params.stressFree = 3; % Width of area in elements left of right side where stress is ignored for L-beam
params.mma = [0.1, 10, 0.01]; % initial values [move, lower, upper] 
params.mmaEnd = [350, 0.05, 0.1, 0.001]; % values after iter [iter, move, lower, upper] 

params.print = [0,0,0];
params.plots = 0;

params.saveName = "";
data = load('input.mat');
fn = fieldnames(data.params);
for k = 1:length(fn) 
    params.(fn{k}) = data.params.(fn{k});
end

sol = Solver(params);
x = ones(sol.nel, 1);
if data.x == 0
    x = 0.5 + (1 - 0.5) * rand(sol.nel, 1);
elseif ~isempty(data.x)
    x = data.x*ones(sol.nel, 1);
end
    
[sol, x] = opt(sol, x);
saveData(sol, x, params);
