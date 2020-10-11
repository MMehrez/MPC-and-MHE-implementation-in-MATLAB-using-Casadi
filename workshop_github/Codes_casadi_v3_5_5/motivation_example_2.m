clear all
close all
clc

x = [0:0.001:4.2*pi];
phi = exp(0.2*x).*sin(x);

Dphi_dt = 0.2*exp(0.2*x).*(sin(x)+5*cos(x));


line_width = 1.5;    fontsize_labels = 15;
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fontsize_labels)

figure
plot(x,phi,'b', 'linewidth',line_width); hold on
plot(x,Dphi_dt,'--k', 'linewidth',line_width); hold on
% plot the minimum points
plot(0,0,'or', 'linewidth',line_width,'MarkerFaceColor','r')
plot(4.9,-2.62,'or', 'linewidth',line_width,'MarkerFaceColor','r')
plot(11.2,-9.2,'or', 'linewidth',line_width,'MarkerFaceColor','r')
legend('\phi(w)','d\phi(w)/dw')
xlabel('w')
ylabel('\phi(w)')
axis([0 14 -10 10])
grid on


%----------------------------------------------
% compose the optimization problem
%----------------------------------------------
% CasADi v3.4.5
% addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.4.5')
% CasADi v3.5.5
addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
 
x = SX.sym('w'); % Decision variables 
obj = exp(0.2*x).*sin(x); % calculate obj

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
args.lbx = 0;  % constrained  
args.ubx = 4*pi;   % constrained 
args.lbg = -inf;  % unconstrained 
args.ubg = inf;   % unconstrained 

args.p   =  [];  % There are no parameters in this optimization problem
args.x0  = 10; % initialization of the optimization problem

sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
    'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
x_sol = full(sol.x)            % Get the solution
min_value = full(sol.f)   % Get the value function





