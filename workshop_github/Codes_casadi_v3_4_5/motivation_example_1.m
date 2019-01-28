clear all
close all
clc

x = [1:0.1:5];
phi = x.^2-6*x+13;


line_width = 1.5;    fontsize_labels = 15;
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fontsize_labels)

figure
plot(x,phi,'b', 'linewidth',line_width); hold on
plot(3,4,'or', 'linewidth',line_width,'MarkerFaceColor','r')
legend('\phi(w)')
xlabel('w')
ylabel('\phi(w)')
axis([0 6 2 10])
grid on


%----------------------------------------------
% compose the optimization problem
%----------------------------------------------
% CasADi v3.4.5
addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.4.5')
import casadi.*
 
x = SX.sym('w'); % Decision variables (controls)
obj = x^2-6*x+13 ; % calculate obj
 
g = [];  % Optimization constraints – empty (unconstrained)
P = [];  % Optimization problem parameters – empty (no parameters used here)
 
OPT_variables = x;  %single decision variable
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);


opts = struct;
opts.ipopt.max_iter = 1000;
opts.ipopt.print_level = 0; %0,3
opts.print_time = 0; %0,1
opts.ipopt.acceptable_tol =1e-8; % optimality convergence tolerance
opts.ipopt.acceptable_obj_change_tol = 1e-6; 

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
args.lbx = -inf;  % unconstrained optimization 
args.ubx = inf;   % unconstrained optimization
args.lbg = -inf;  % unconstrained optimization
args.ubg = inf;   % unconstrained optimization

args.p   =  [];  % There are no parameters in this optimization problem
args.x0  = -0.5; % initialization of the optimization problem

sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
    'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
x_sol = full(sol.x)            % Get the solution
min_value = full(sol.f)   % Get the value function





