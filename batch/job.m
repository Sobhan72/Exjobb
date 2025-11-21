clc, clear, close all
%% Add path to data
addpath(genpath('/home/zachariand/Exjobb'));
addpath(genpath('/home/zachariand/Exjobb/functions'));
%% Run Job
% FEM parameters

params.le = 0.002;
params.lx = 0.16; 
params.ly = 0.16; 
params.wx = [];
params.wy = []; 
params.loadcase = 4;

params.t = 0.1;
params.ngp = 4;

params.R1tol = 1e-2;
% params.N = 4;
params.disp = -1.6e-3; % Total displacement 

% Material parameters
E = 210e9;
v = 0.3;
params.E1 = E; params.E2 = 0.7*E; params.E3 = 0.5*E;
params.v12 = v; params.v13 = 0.5*v; params.v23 = 0.5*v;

params.sigy01 = 360e6;
params.sigy02 = 300e6;
params.sigy03 = 250e6;

params.H = 10e9;
params.Kinf = 0; %0.3*params.sigy01; %extra terms linHard (sat. stress)
params.xi = 1e-3; %extra terms linHard (sat. exp)

params.rtol = 1e-4;
params.PT = 1; % 1 for Deformation plasticity, 0 for Elastoplasticity

params.ngr = 1;

% Optimization Parameters
params.re = 3; % Elements in radius (dubblerad)
params.filtOn = true;
params.pad = true;

params.p = 3;
params.q = 2.5;
params.del = 1e-9;
params.dels = 1e-3;
params.rampPQ = [4, 0.1+2/30]; % [end value of p, increment size]

params.eta = 0.5;
params.beta = 0.1;
params.rampB = [2, 10, 1.13]; % [0/1/2, end value, factor size]  (0: off, 1: on, 2: on after simp converges)

%params.ncon = 1; % Nr of constraints
params.Vf = 0.3;
params.xtol = 1e-3;
params.ftol = 0.2;
% params.iterMax = 1000;

params.stressCon = 0;
params.pnm = 8;
params.sigc = 1.05; %360e6; % Max stress for constraint
params.stressFree = 0; % Width of area in elements left of right side where stress is ignored for L-beam
params.zeroGrad = false; % Manually zero stress gradient in stressFree zone.
params.plasticFree = 0; % Width of area in elements left of right boundary where plasticity is ignored for L-beam
params.mma = [0.1, 10, 0.01]; % initial values [move, lower, upper] 
params.mmaEnd = [0, 0.05, 0.1, 0.001]; % values after iter 350 [iter, move, lower, upper] 

params.ang = 45;

params.prints = [0,0,0]; %[0,0,0];
params.plots = 0;

params.saveName = "";
data = load('input1.mat');
fn = fieldnames(data.params);
for k = 1:length(fn) 
    params.(fn{k}) = data.params.(fn{k});
end
params.N = 3;
params.iterMax = 450;
sol = Solver(params);
x = ones(sol.nel, 1);
if data.x == 0
    x = 0.5 + (1 - 0.5) * rand(sol.nel, 1);
elseif ~isempty(data.x)
    x = data.x*ones(sol.nel, 1);
end
    
[sol, x] = opt(sol, x);
saveData(sol, x, params);
