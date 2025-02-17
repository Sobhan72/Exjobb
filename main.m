clc, clear, close all

% Input parameters
params.le = 1;
params.lx = 4;
params.ly = 4;

params.E = 210e9;
params.v = 0.3;
params.ep = [2 1 2];
params.sig_y0 = 360e6;
params.H = 10e9;

params.Fco = 1/(2*params.sig_y0^2); 
params.Gco = 1/(2*params.sig_y0^2); 
params.Hco = 1/(2*params.sig_y0^2); 
params.Lco = 3/(2*params.sig_y0^2);


sol = Solver(params);
patch(sol.ex', sol.ey', 1)