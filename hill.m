%Hill elastoplastic
clc, clear all, close all
E = 210e9; v = 0.3; sig_y0 = 360e6; Hard = 10e9; K = E/(3*(1-2*v));
F = 1/(2*sig_y0); G = 1/(2*sig_y0); H = 1/(2*sig_y0);
L = 3/(2*sig_y0); M = 3/(2*sig_y0); N = 3/(2*sig_y0);
P = [F+G -F -G 0 0 0; -F F+H -H 0 0 0; -G -H G+H 0 0 0; 0 0 0 2*L 0 0; 0 0 0 0 2*M 0; 0 0 0 0 0 2*N];