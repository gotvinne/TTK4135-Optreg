% TTK4135 - Helicopter lab
% Updated spring 2018, Andreas L. Flåten

%Script calculating optimal trajectory for the given system
%In addition to model, we optimise with inequality constraint. 
%x = [lambda, r, p, p_dot, e, e_dot]'

init; % Change this to the init file corresponding to your helicopter

% Initialization and model definition

% Discrete time system model. x = [lambda r p p_dot]'
T = 0.25; % sampling time
Ad = [1,T,0,0,0,0;
    0,1,-T*K_2,0,0,0;
    0,0,1,T,0,0;
    0,0,-T*K_1*K_pp,1-T*K_1*K_pd,0,0;
    0,0,0,0,1,T;
    0,0,0,0,-T*K_3*K_ep,1-T*K_3*K_ed];

Bd = [0,0;
    0,0;
    0,0;
    T*K_1*K_pp,0;
    0,0;
    0,T*K_3*K_ep];

% Number of states and inputs
mx = size(Ad,2); % Number of states (number of columns in A)
mu = size(Bd,2); % Number of inputs(number of columns in B)

%% LQR
R_elements = [0.1, 0.01];
Q_elements = [5, 1, 1, 1, 10, 1];

R  = diag(R_elements);
Q = diag(Q_elements);

[K, S, e] = dlqr(Ad,Bd,Q,R);


%% Optimalization problem 

% Initial values
lambda_0 = pi;
lambda_f = 0;

x1_0 = lambda_0;                        % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               %e
x6_0 = 0;                               %e_dot
x0 = [x1_0 x2_0 x3_0 x4_0,x5_0,x6_0]';           % Initial values

% Time horizon and initialization
N  = 60;                                 % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon,
z0 = z;             % Initial value for optimization
z0(1:6) = x0;

% Bounds
ul 	    = -(30*pi)/180;                   % Lower bound on control
uu 	    = (30*pi)/180;                   % Upper bound on control

xl      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul;                           % Lower bound on state x3, korresponding to pc.
xu(3)   = uu;                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu); % hint: gen_constraints
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
Q1(5,5) = 0;                            % Weight on state x5
Q1(6,6) = 0;                            % Weight on state x6

q1 = 1; 
q2 = 1; 
R = diag([q1,q2]);                            
G = gen_q(Q1,R,N,M);                       % Generate Q, hint: gen_q

% Generate system matrixes for linear model
Aeq = gen_aeq(Ad,Bd,N,mx,mu);             % Generate A, hint: gen_aeq
beq = zeros(size(Aeq,1),1); % Generate b
beq(1:6) = Ad*x0;


%% Solving cost function: 

B = 0; 
options = optimset('Algorithm','sqp','MaxFunEvals',40000); 

fun = @(z)(1/2)*z.'*G*z; 
z = fmincon(fun,z0,[],[],Aeq,beq,vlb,vub,@inequality_constraint, options);


%% Extract optimalization variable: 

u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution
pc = [u(1:mu:N*mu);u(M*mu-1)]; 
ec = [u(2:mu:N*mu);u(M*mu)];

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution

num_variables = 5/T;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);


%Optimal solution for horizon N:
pc   = [zero_padding; pc];
ec   = [zero_padding; ec];
x1  = [pi*unit_padding; x1];
x2  = [zero_padding; x2];
x3  = [zero_padding; x3];
x4  = [zero_padding; x4];
x5 = [zero_padding; x5]; 
x6 = [zero_padding; x6]; 
time_steps = [0:T:(size(x1, 1) - 1) * T]';

 
save('simulated_data.mat','time_steps','pc','ec','x1','x2','x3','x4','x5','x6');

pc = timeseries(pc,time_steps); 
ec = timeseries(ec,time_steps);
x = [x1'; x2'; x3'; x4'; x5'; x6'];
x = timeseries(x, time_steps);

% optimal_trajectory_experimental_plot;

function [c,ceq] = inequality_constraint(z)
    N = 60;
    a = 0.2; 
    b = 20;
    
    mx = 6; 
    lambda_t = 2*pi/3;
    c = zeros(N,1);
    for i = 1:N
        c(i) = a*exp(-b*(z((i-1)*mx+1)-lambda_t)^2)-z((i-1)*mx+5); 
    end 
    ceq = []; 
end

