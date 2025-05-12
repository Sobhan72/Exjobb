clc, clear, close all

load("input.mat");

% FEM parameters
params.le = 0.01;
params.lx = 0.1;
params.ly = 0.05;
params.Vf = 0.3;

params.t = 1;
params.ngp = 4;

params.R1tol = 1e-1;
params.N = 4;
% params.disp = -7e-4; % displacement [nodes total-size]

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
params.p = 1.5;
params.q = 1;
params.eta = 0.5;
params.beta = 1;
params.ramp = true;
params.del = 1e-9;
params.ncon = 1; % Nr of constraints
params.xtol = 1e-5;
params.iterMax = 10;


params.saveName = "OptDesign1";
sol = Solver(params);

x = ones(sol.nel, 1);
[sol, x] = optimizer(sol, x);
saveData(sol, x, params);
