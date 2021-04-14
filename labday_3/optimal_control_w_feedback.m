% TTK4135 - Helicopter lab
% Updated spring 2018, Andreas L. Flåten

%Script solving task given in lab 3. 

init; % Change this to the init file corresponding to your helicopter

% Initialization and model definition


% Discrete time system model. x = [lambda r p p_dot]'
T = 0.25; % sampling time
Ad = [1, T, 0, 0;
    0, 1, -T*K_2, 0;
    0, 0, 1, T;
    0, 0, -T*K_1*K_pd, 1- T*K_1*K_pd];

Bd = [0;
    0;
    0;
    T*K_1*K_pp];

% Number of states and inputs
mx = size(Ad,2); % Number of states (number of columns in A)
mu = size(Bd,2); % Number of inputs(number of columns in B)

%% LQR matrices

R_elements = [1];
Q_elements = [1,1,1,1];

R  = diag(R_elements);
Q = diag(Q_elements);

%Calculate optimal gain matrix K
%K=(BTSB+R)−1(BTSA+NT)
[K,S,e] = dlqr(Ad,Bd,Q,R);


% Initial values
lambda_0 = pi;
lambda_f = 0;

x1_0 = lambda_0;                               % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x0 = [x1_0 x2_0 x3_0 x4_0]';           % Initial values

% Time horizon and initialization
N  = 100;                                 % Time horizon for states
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon,
z0 = z;             % Initial value for optimization
z0(1:4) = x0;

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

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem)
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4

P1 = 1;                                % Weight on input, try 0.1, 1 and 10
Q = gen_q(Q1,P1,N,M);                                  % Generate Q, hint: gen_q
c = 0;                                  % Generate c, this is the linear constant term in the QP, due to no linear terms

% Generate system matrixes for linear model
Aeq = gen_aeq(Ad,Bd,N,mx,mu);             % Generate A, hint: gen_aeq
beq = zeros(size(Aeq,1),1); % Generate b
beq(1:4) = Ad*x0;


% Solve QP problem with linear model
tic

%zlb = [repmat(xl,N,1); repmat(ul,N,1)];
%zub = [repmat(xu,N,1); repmat(uu,N,1)];


%% Solving cost function: 

[z,fval] = quadprog(Q,[],[],[],Aeq,beq,vlb,vub); % hint: quadprog. Type 'doc quadprog' for more info
t1=toc;

% Calculate objective value
phi1 = 0.0;
PhiOut = zeros(N*mx+M*mu,1);
for i=1:N*mx+M*mu
    phi1=phi1+Q(i,i)*z(i)*z(i);
    PhiOut(i) = phi1;
end

% Extract control inputs and states
u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution

num_variables = 5/T;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

%Optimal solution for horizon N:
u   = [zero_padding; u; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];

%% To simulink:

% Added 5 seconds of zeros on both sides of the data
time_steps = [0:T:N * T + 5 * 2]';

u = timeseries(u,time_steps); 
x_state = [x1';x2';x3';x4']; 
x = timeseries(x_state,time_steps); 

