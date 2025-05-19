clc, clear, close all
%% Add path to data
addpath(genpath('/home/zachariand/Exjobb'));
addpath(genpath('/home/zachariand/Exjobb/functions'));
%% Run Job
% FEM parameters
params.le = 0.001;
params.lx = 0.1;
params.ly = 0.05;
params.Vf = 0.3;

params.t = 1;
params.ngp = 4;

params.R1tol = 1e-1;
params.N = 3;
params.disp = -1e-3; % displacement [nodes total-size]

% Material parameters
E = 210e9;
v = 0.3;
params.E1 = E; params.E2 = E; params.E3 = E;
params.v12 = v; params.v13 = v; params.v23 = v;
params.v21 = v; params.v32 = v; params.v31 = v;

params.sigy01 = 360e6;
params.sigy02 = 360e6;
params.sigy03 = 360e6;
 
params.H = 10e9;
params.Kinf = 0; %0.3*params.sigy01; %extra terms linHard (sat. stress)
params.xi = 1e-3; %extra terms linHard (sat. exp)

params.rtol = 1e-1;
params.DP = 1; % 1 for Deformation plasticity, 0 for Elastoplasticity

% Optimization Parameters
params.re = 3; % Elements in radius
params.filtOn = true;
params.loadcase = 1;
params.p = 1.5; %3;
params.q = 1; %2,5;
params.eta = 0.5;
params.beta = 1;
params.rampB = 1; % 0: off, 1: B*1.1, 2: B + 1
params.rampPQ = true;
params.del = 1e-9;
params.ncon = 1; % Nr of constraints
params.xtol = 1e-5;
params.iterMax = 500;

params.saveName = "";
data = load('input.mat');
fn = fieldnames(data.params);
for k = 1:length(fn) 
    params.(fn{k}) = data.params.(fn{k});
end

sol = Solver(params);
x = ones(sol.nel, 1);
if ~isempty(data.x)
    x = data.x*ones(sol.nel, 1);
end
[sol, x] = optimizer(sol, x);
saveData(sol, x, params);

