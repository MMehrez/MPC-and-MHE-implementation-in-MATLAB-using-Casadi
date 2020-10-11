% first casadi test for mpc fpr mobile robots
clear all
close all
clc

% CasADi v3.1.1
% addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.4.5')
% CasADi v3.5.5
addpath('C:\Users\mehre\OneDrive\Desktop\CasADi\casadi-windows-matlabR2016a-v3.5.5')
import casadi.*

T = 0.2; %[s]
N = 10; % prediction horizon
rob_diam = 0.3;

v_max = 0.6; v_min = -v_max;
omega_max = pi/4; omega_min = -omega_max;

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);

v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s

f = Function('f',{states,controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U',n_controls,N); % Decision variables (controls)
P = SX.sym('P',n_states + n_states);
% parameters (which include at the initial state of the robot and the reference state)

X = SX.sym('X',n_states,(N+1));
% A vector that represents the states over the optimization problem.

obj = 0; % Objective function
g = [];  % constraints vector

Q = zeros(3,3); Q(1,1) = 1;Q(2,2) = 5;Q(3,3) = 0.1; % weighing matrices (states)
R = zeros(2,2); R(1,1) = 0.5; R(2,2) = 0.05; % weighing matrices (controls)

st  = X(:,1); % initial state
g = [g;st-P(1:3)]; % initial condition constraints
for k = 1:N
    st = X(:,k);  con = U(:,k);
    obj = obj+(st-P(4:6))'*Q*(st-P(4:6)) + con'*R*con; % calculate obj
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end
% make the decision variable one column  vector
OPT_variables = [reshape(X,3*(N+1),1);reshape(U,2*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);


args = struct;

args.lbg(1:3*(N+1)) = 0;
args.ubg(1:3*(N+1)) = 0;

args.lbx(1:3:3*(N+1),1) = -2; %state x lower bound
args.ubx(1:3:3*(N+1),1) = 2; %state x upper bound
args.lbx(2:3:3*(N+1),1) = -2; %state y lower bound
args.ubx(2:3:3*(N+1),1) = 2; %state y upper bound
args.lbx(3:3:3*(N+1),1) = -inf; %state theta lower bound
args.ubx(3:3:3*(N+1),1) = inf; %state theta upper bound

args.lbx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_min; %v lower bound
args.ubx(3*(N+1)+1:2:3*(N+1)+2*N,1) = v_max; %v upper bound
args.lbx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_min; %omega lower bound
args.ubx(3*(N+1)+2:2:3*(N+1)+2*N,1) = omega_max; %omega upper bound
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0.1 ; 0.1 ; 0.0];    % initial condition.
xs = [1.5 ; 1.5 ; 0.0]; % Reference posture.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,2);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 20; % total sampling times

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
tic
while(norm((x0-xs),2) > 0.05 && mpciter < sim_tim / T)
    args.p   = [x0;xs]; % set the values of the parameters vector
    % initial value of the optimization variables
    args.x0  = [reshape(X0',3*(N+1),1);reshape(u0',2*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x(3*(N+1)+1:end))',2,N)'; % get controls only from the solution
    xx1(:,1:3,mpciter+1)= reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    % Apply the control and shift the solution
    [t0, x0, u0] = shift(T, t0, x0, u,f);
    xx(:,mpciter+2) = x0;
    X0 = reshape(full(sol.x(1:3*(N+1)))',3,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:);X0(end,:)];
    mpciter
    mpciter = mpciter + 1;
end;
toc

ss_error = norm((x0-xs),2)
Draw_MPC_point_stabilization_v1 (t,xx,xx1,u_cl,xs,N,rob_diam)
%-----------------------------------------
%-----------------------------------------
%-----------------------------------------
%    Start MHE implementation from here
%-----------------------------------------
%-----------------------------------------
%-----------------------------------------
% plot the ground truth
figure(1)
subplot(311)
plot(t,xx(1,1:end-1),'b','linewidth',1.5); axis([0 t(end) 0 1.8]);hold on
ylabel('x (m)')
grid on
subplot(312)
plot(t,xx(2,1:end-1),'b','linewidth',1.5); axis([0 t(end) 0 1.8]);hold on
ylabel('y (m)')
grid on
subplot(313)
plot(t,xx(3,1:end-1),'b','linewidth',1.5); axis([0 t(end) -pi/4 pi/2]);hold on
xlabel('time (seconds)')
ylabel('\theta (rad)')
grid on

% Synthesize the measurments
con_cov = diag([0.005 deg2rad(2)]).^2;
meas_cov = diag([0.1 deg2rad(2)]).^2;

r = [];
alpha = [];
for k = 1: length(xx(1,:))-1
    r = [r; sqrt(xx(1,k)^2+xx(2,k)^2)  + sqrt(meas_cov(1,1))*randn(1)];
    alpha = [alpha; atan(xx(2,k)/xx(1,k))      + sqrt(meas_cov(2,2))*randn(1)];
end
y_measurements = [ r , alpha ];

% Plot the cartesian coordinates from the measurements used
figure(1)
subplot(311)
plot(t,r.*cos(alpha),'r','linewidth',1.5); hold on
grid on
legend('Ground Truth','Measurement')
subplot(312)
plot(t,r.*sin(alpha),'r','linewidth',1.5); hold on
grid on

% plot the ground truth mesurements VS the noisy measurements
figure(2)
subplot(211)
plot(t,sqrt(xx(1,1:end-1).^2+xx(2,1:end-1).^2),'b','linewidth',1.5); hold on
plot(t,r,'r','linewidth',1.5); axis([0 t(end) -0.2 3]); hold on
ylabel('Range: [ r (m) ]')
grid on
legend('Ground Truth','Measurement')
subplot(212)
plot(t,atan(xx(2,1:end-1)./xx(1,1:end-1)),'b','linewidth',1.5); hold on
plot(t,alpha,'r','linewidth',1.5); axis([0 t(end) 0.2 1]); hold on
ylabel('Bearing: [ \alpha (rad) ]')
grid on

% The following two matrices contain what we know about the system, i.e.
% the nominal control actions applied (measured) and the range and bearing
% measurements.
%-------------------------------------------------------------------------
u_cl;
y_measurements;

% Let's now formulate our MHE problem
%------------------------------------
T = 0.2; %[s]
N_MHE = size(y_measurements,1)-1; % Estimation horizon

v_max = 0.6; v_min = -v_max;
omega_max = pi/4; omega_min = -omega_max;

x = SX.sym('x'); y = SX.sym('y'); theta = SX.sym('theta');
states = [x;y;theta]; n_states = length(states);
v = SX.sym('v'); omega = SX.sym('omega');
controls = [v;omega]; n_controls = length(controls);
rhs = [v*cos(theta);v*sin(theta);omega]; % system r.h.s
f = Function('f',{states,controls},{rhs}); % MOTION MODEL

r = SX.sym('r'); alpha = SX.sym('alpha'); % range and bearing
measurement_rhs = [sqrt(x^2+y^2); atan(y/x)];
h = Function('h',{states},{measurement_rhs}); % MEASUREMENT MODEL
%y_tilde = [r;alpha];

% Decision variables
U = SX.sym('U',n_controls,N_MHE);    %(controls)
X = SX.sym('X',n_states,(N_MHE+1));  %(states) [remember multiple shooting]

P = SX.sym('P', 2 , N_MHE + (N_MHE+1));
% parameters (include r and alpha measurements as well as controls measurements)

V = inv(sqrt(meas_cov)); % weighing matrices (output)  y_tilde - y
W = inv(sqrt(con_cov)); % weighing matrices (input)   u_tilde - u

obj = 0; % Objective function
g = [];  % constraints vector
for k = 1:N_MHE+1
    st = X(:,k);
    h_x = h(st);
    y_tilde = P(:,k);
    obj = obj+ (y_tilde-h_x)' * V * (y_tilde-h_x); % calculate obj
end

for k = 1:N_MHE
    con = U(:,k);
    u_tilde = P(:,N_MHE+ k);
    obj = obj+ (u_tilde-con)' * W * (u_tilde-con); % calculate obj
end

% multiple shooting constraints
for k = 1:N_MHE
    st = X(:,k);  con = U(:,k);
    st_next = X(:,k+1);
    f_value = f(st,con);
    st_next_euler = st+ (T*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end


% make the decision variable one column  vector
OPT_variables = [reshape(X,3*(N_MHE+1),1);reshape(U,2*N_MHE,1)];

nlp_mhe = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_mhe,opts);

args = struct;

args.lbg(1:3*(N_MHE)) = 0; % equality constraints
args.ubg(1:3*(N_MHE)) = 0; % equality constraints

args.lbx(1:3:3*(N_MHE+1),1) = -2; %state x lower bound
args.ubx(1:3:3*(N_MHE+1),1) = 2; %state x upper bound
args.lbx(2:3:3*(N_MHE+1),1) = -2; %state y lower bound
args.ubx(2:3:3*(N_MHE+1),1) = 2; %state y upper bound
args.lbx(3:3:3*(N_MHE+1),1) = -pi/2; %state theta lower bound
args.ubx(3:3:3*(N_MHE+1),1) = pi/2; %state theta upper bound

args.lbx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE,1) = v_min; %v lower bound
args.ubx(3*(N_MHE+1)+1:2:3*(N_MHE+1)+2*N_MHE,1) = v_max; %v upper bound
args.lbx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE,1) = omega_min; %omega lower bound
args.ubx(3*(N_MHE+1)+2:2:3*(N_MHE+1)+2*N_MHE,1) = omega_max; %omega upper bound
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP

% MHE Simulation
%------------------------------------------
U0 = zeros(N_MHE,2);   % two control inputs for each robot
X0 = zeros(N_MHE+1,3); % initialization of the states decision variables

U0 = u_cl(1:N_MHE,:); % initialize the control actions by the measured
% initialize the states from the measured range and bearing
X0(:,1:2) = [y_measurements(1:N_MHE+1,1).*cos(y_measurements(1:N_MHE+1,2)),...
    y_measurements(1:N_MHE+1,1).*sin(y_measurements(1:N_MHE+1,2))];

k=1;
% Get the measurements window and put it as parameters in p
args.p   = [y_measurements(k:k+N_MHE,:)',u_cl(k:k+N_MHE-1,:)'];
% initial value of the optimization variables
args.x0  = [reshape(X0',3*(N_MHE+1),1);reshape(U0',2*N_MHE,1)];
sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
    'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
U_sol = reshape(full(sol.x(3*(N_MHE+1)+1:end))',2,N_MHE)'; % get controls only from the solution
X_sol = reshape(full(sol.x(1:3*(N_MHE+1)))',3,N_MHE+1)'; % get solution TRAJECTORY

figure(1)
subplot(311)
plot(t,X_sol(:,1),'g','linewidth',1.5); hold on
legend('Ground Truth','Measurement','MHE')
subplot(312)
plot(t,X_sol(:,2),'g','linewidth',1.5); hold on
subplot(313)
plot(t,X_sol(:,3),'g','linewidth',1.5); hold on


