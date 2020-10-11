clear all
close all
clc

x = [0,45,90,135,180];
y = [667,661,757,871,1210];

line_width = 1.5;    fontsize_labels = 15;
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', fontsize_labels)

figure(1)
plot(x,y,'*b', 'linewidth',line_width); hold on
xlabel('x')
ylabel('y')
grid on


%----------------------------------------------
% compose the optimization problem
%----------------------------------------------
% CasADi v3.4.5
% addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.4.5')
% CasADi v3.5.5
addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*
 
m = SX.sym('m'); % Decision variable (slope)
c = SX.sym('c'); % Decision variable (y-intersection)

obj = 0;

for i = 1:length(x)
   obj = obj+ (y(i) - (m*x(i)+c))^2; 
end

g = [];   % Optimization constraints – empty (unconstrained)
P = [];  % Optimization problem parameters – empty (no parameters used here)
 
OPT_variables = [m,c];  %Two decision variable
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
args.x0  = [0.5,1]; % initialization of the optimization problem

sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
    'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
x_sol = full(sol.x);            % Get the solution
min_value = full(sol.f)   % Get the value function

x_line = [0:1:180];
m_sol = x_sol(1)
c_sol = x_sol(2)
y_line = m_sol*x_line+c_sol;
figure(1)
plot(x_line,y_line,'-k', 'linewidth',line_width); hold on
legend('Data points','y = 2.88 \cdot x + 574')


% visualize the objective function
obj_fun = Function('obj_fun',{m,c},{obj});

m_range = [-1:0.5:6];
c_range = [400:50:800];
obj_plot_data = [];

[mm,cc] = meshgrid(m_range,c_range);

for n = 1:1:size(mm,1) 
    for k = 1:1:size(mm,2) %
      obj_plot_data(n,k) = full(obj_fun(mm(n,k),cc(n,k))) ;
    end
end

figure
surf(mm,cc,obj_plot_data); hold on
xlabel('(m)')
ylabel('(c)')
zlabel('(\phi)')

box on
ax = gca;
ax.BoxStyle = 'full';

min(min(obj_plot_data))
