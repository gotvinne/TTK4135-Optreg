close all;

data = load('data.mat').ans;
t_s = load('simulated_data.mat').time_steps;
p_c_s = load('simulated_data.mat').pc' * 180/pi;
e_c_s = load('simulated_data.mat').ec' * 180/pi;
travel_s = load('simulated_data.mat').x1 * 180/pi;
travel_rate_s = load('simulated_data.mat').x2 * 180/pi;
pitch_s = load('simulated_data.mat').x3 * 180/pi;
pitch_rate_s = load('simulated_data.mat').x4 * 180/pi;
elevation_s = load('simulated_data.mat').x5 * 180/pi;
elevation_rate_s = load('simulated_data.mat').x6 * 180/pi;
T = 0.25;

t = data(1, :);
p_c = data(2, :) * 180 / pi;
e_c = data(3, :) * 180 / pi;
travel = data(4, :);
travel_rate = data(5, :);
pitch = data(6, :);
pitch_rate = data(7, :);
elevation = data(8, :);
elevation_rate = data(9, :);

%% Plotting 2.0 
    

    % Added 5 seconds of zeros on both sides of the data

    figure;
    
    subplot(421);
    title('Optimal trajectory');
    plot(t,travel);
    hold on;
    plot(t_s,travel_s);
    grid on;
    legend('Experimental', 'Simulated');
    ylabel('lambda')
    
    subplot(422);
    plot(t,travel_rate);
    hold on;
    plot(t_s,travel_rate_s);
    grid on;
    legend('Experimental', 'Simulated');
    ylabel('lambda_{dot} [rad/s]')
    
    subplot(423) 
    plot(t,pitch);
    hold on;
    plot(t_s,pitch_s);
    grid on;
    legend('Experimental', 'Simulated');
    ylabel('p [rad]')
    
    subplot(424)
    plot(t,pitch_rate);
    hold on;
    plot(t_s,pitch_rate_s);
    grid on;
    legend('Experimental', 'Simulated');
    ylabel('p_{dot} [rad/s]')
    
    
    subplot(425)
    plot(t,elevation);
    hold on;
    plot(t_s,elevation_s);
    hold on;
    
    constraint_values = [];
    for lambda = travel_s'
        constraint_values = [constraint_values constraint(lambda * pi / 180) * 180 / pi];
    end
    
    plot(t_s, constraint_values);
    grid on;
    
    legend('Experimental', 'Simulated', 'Constraint');
    xlabel('time (s)'),ylabel('e [rad]')
    
    subplot(426)
    plot(t,elevation_rate);
    hold on;
    plot(t_s,elevation_rate_s);
    grid on;
    legend('Experimental', 'Simulated');
    xlabel('time (s)'),ylabel('e_{dot} [rad/s]')
    
    subplot(427)
    plot(t,p_c);
    hold on;
    plot(t_s,p_c_s);
    grid on;
    legend('Experimental', 'Simulated');
    xlabel('time (s)'),ylabel('p_c [rad]')
    
    subplot(428)
    plot(t,e_c);
    hold on;
    plot(t_s,e_c_s);
    grid on; 
    legend('Experimental', 'Simulated');
    xlabel('time (s)'),ylabel('e_c [rad]')

    print('data.eps', '-depsc');
    
function elevation = constraint(lambda)
    a = 0.2; 
    b = 20;
    lambda_t = 2*pi/3;
    elevation = a*exp(-b*(lambda-lambda_t)^2);
end