
data = load('data.mat').ans;
t_s = load('simulated_data.mat').time_steps;
u_s = load('simulated_data.mat').u';
travel_s = load('simulated_data.mat').x1;
travel_rate_s = load('simulated_data.mat').x2;
pitch_s = load('simulated_data.mat').x3;
pitch_rate_s = load('simulated_data.mat').x4;
T = 0.25;
padding_time = 5;

padding = zeros(padding_time / T, 1);

travel_s = [ones(padding_time / T, 1) * pi; travel_s; padding];
travel_rate_s = [padding; travel_rate_s; padding];
pitch_s = [padding; pitch_s; padding];
pitch_rate_s = [padding; pitch_rate_s; padding];

t = data(1, :);
u = data(2, :);
travel = data(3, :);
travel_rate = data(4, :);
pitch = data(5, :);
pitch_rate = data(6, :);


x1_array = travel * pi/180 + pi;
x2_array = travel_rate * pi/180;
x3_array = pitch * pi/180;
x4_array = pitch_rate * pi/180;
u_array = u; 
t_plot = 0; 
xEnd = 17;

%% Plotting 2.0 
    
    figure(3)
    
    subplot(511)
    title('Optimal trajectory, experimental');
    hold on;
    ylim([-1 1]);
    xlim([0 xEnd]);
    plot(t,u_array),grid
    plot(t_s,u_s),grid
    hold off; 
    legend('Experimental','Simulated');
    ylabel('u [rad]')
    
    subplot(512)
    
    hold on; 
    xlim([0 xEnd]);
    plot(t,x1_array),grid
    plot(t_s,travel_s),grid
    hold off; 
    legend('Experimental','Simulated');
    ylabel('lambda [rad]')
    
    subplot(513)
    hold on; 
    xlim([0 xEnd]);
    plot(t,x2_array),grid
    plot(t_s,travel_rate_s),grid
    hold off; 
    legend('Experimental','Simulated');
    ylabel('r [rad/s]')
    
    subplot(514)
    hold on; 
    xlim([0 xEnd]);
    plot(t,x3_array),grid
    plot(t_s,pitch_s),grid
    ylabel('p [rad]')
    hold off; 
    legend('Experimental','Simulated');
    
    subplot(515)
    hold on;
    xlim([0 xEnd]);
    plot(t,x4_array),grid
    plot(t_s,pitch_rate_s),grid
    hold off; 
    legend('Experimental','Simulated');
    xlabel('time (s)'),ylabel('pdot [rad/s]')
    
    print('data.eps', '-depsc');
    