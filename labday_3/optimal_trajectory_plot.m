% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Flåten

init; % Change this to the init file corresponding to your helicopter
q = [0.1, 1, 10]; 
x1_array = [];
x2_array = [];
x3_array = [];
x4_array = [];
u_array = []; 
t_plot = 0; 

iteration = 1; 
for iteration = 1:length(q)  

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

    
    P1 = q(iteration);                                % Weight on input, try 0.1, 1 and 10 
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


    [z,fval] = quadprog(Q,[],[],[],Aeq,beq,vlb,vub); % hint: quadprog. Type 'doc quadprog' for more info 
    t1=toc;

    % Calculate objective value
    %{phi1 = 0.0;
    
%     PhiOut = zeros(N*mx+M*mu,1);
%     for i=1:N*mx+M*mu
%         phi1=phi1+Q(i,i)*z(i)*z(i);
%         PhiOut(i) = phi1;
%     end
    %}

    % Extract control inputs and states
    u  = [z(N*mx+1:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution

    x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
    x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
    x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
    x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution

    num_variables = 5/T;
    zero_padding = zeros(num_variables,1);
    unit_padding  = ones(num_variables,1);

    u   = [zero_padding; u; zero_padding];
    x1  = [pi*unit_padding; x1; zero_padding];
    x2  = [zero_padding; x2; zero_padding];
    x3  = [zero_padding; x3; zero_padding];
    x4  = [zero_padding; x4; zero_padding];
    
    x1_array = [x1_array, x1];
    x2_array = [x2_array, x2]; 
    x3_array = [x3_array, x3]; 
    x4_array = [x4_array, x4]; 
    u_array = [u_array, u];

    % Plotting
    
    t = 0:T:T*(length(u)-1);
    %{
    color = ['r','b','g']; 

    figure(2)
    subplot(511)
    stairs(t,u),grid
    ylabel('u')
    subplot(512)
    plot(t,x1,'m',t,x1,'mo','Color', color(iteration)),grid
    ylabel('lambda')
    subplot(513)
    plot(t,x2,'m',t,x2','mo','Color', color(iteration)),grid
    ylabel('r')
    subplot(514)
    plot(t,x3,'m',t,x3,'mo','Color', color(iteration)),grid
    ylabel('p')
    subplot(515)
    plot(t,x4,'m',t,x4','mo','Color', color(iteration)),grid
    xlabel('tid (s)'),ylabel('pdot')
    %}

end

%% Plotting 2.0 
    
    figure(3)
    
    subplot(511)
    title('Optimal trajectory');
    hold on;
    stairs(t,u_array(:,1)),grid
    stairs(t,u_array(:,2)),grid
    stairs(t,u_array(:,3)),grid
    hold off; 
    legend('0.1','1','10');
    ylabel('u [rad]')
    
    subplot(512)
    
    hold on; 
    plot(t,x1_array(:,1)),grid
    plot(t,x1_array(:,2)),grid
    plot(t,x1_array(:,3)),grid
    hold off; 
    legend('0.1','1','10');
    ylabel('lambda [rad]')
    
    subplot(513)
    hold on; 
    plot(t,x2_array(:,1)),grid
    plot(t,x2_array(:,2)),grid
    plot(t,x2_array(:,3)),grid
    hold off; 
    legend('0.1','1','10');
    ylabel('r [rad/s]')
    
    subplot(514)
    hold on; 
    plot(t,x3_array(:,1)),grid
    plot(t,x3_array(:,2)),grid
    plot(t,x3_array(:,3)),grid
    ylabel('p [rad]')
    hold off; 
    legend('0.1','1','10');
    
    subplot(515)
    hold on;
    plot(t,x4_array(:,1)),grid
    plot(t,x4_array(:,2)),grid
    plot(t,x4_array(:,3)),grid
    hold off; 
    legend('0.1','1','10');
    xlabel('tid (s)'),ylabel('pdot [rad/s]')
    

