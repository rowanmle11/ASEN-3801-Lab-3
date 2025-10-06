% Zeena Ibrahim
% ASEN 3801 Lab 3
% Date 10/04/2025

clc;
clear;
close all;


% Lab 3.2 Code

trial_torque4 = "RWHEEL_T4t5.xlsx";
trial_torque8 = "RWHEEL_T8t5.xlsx";
trial_torque12 = "RWHEEL_T12t5.xlsx";
trial_torque16 = "RWHEEL_T16t5.xlsx";
trial_torque20 = "RWHEEL_T20t5.xlsx";

trials = [trial_torque4 trial_torque8 trial_torque12 trial_torque16 trial_torque20];
trials_data = cell(1, length(trials));
torques = [4 8 12 16 20];

for i = 1:length(trials)

    trials_data{i} = readmatrix(trials(i));
    time{i} = (trials_data{i}(:,1))/(10^3); %conversion from ms to sec
    torque{i} = (trials_data{i}(:,4))*(33.5); %multiply by torque constant
    angular_velocity{i} = (trials_data{i}(:,3))* ((2*pi)/60); %Convert from RPM to rad/s

     if i <= 3
        t_min = 1; t_max = 6;
    elseif i == 4
        t_min = 1; t_max = 5;
    else
        t_min = 1; t_max = 4;
    end

    linear_idx = (time{i} >= t_min & time{i} <= t_max);
    p = polyfit(time{i}(linear_idx), angular_velocity{i}(linear_idx), 1);
    w_fit = polyval(p, time{i}(linear_idx));
   
    figure()
    plot(time{i}, angular_velocity{i}, 'b');
    hold on;
    plot(time{i}(linear_idx), w_fit, 'r', 'LineWidth', 2);
    title(['Angular Velocity for Trial ' num2str(i) ': ' num2str(torques(i)) ' Nm']);
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    legend('Measured Data', 'Linear Fit', 'Location', 'best');

    alpha(i) = p(1);
end

moi = torques./alpha;

stdev = std(moi);
mean_moi = mean(moi); 
fprintf('Standard Deviation of Moment of Inertia: %.2f kg*m^2\n', stdev);
fprintf('Mean Moment of Inertia: %.2f kg*m^2\n', mean_moi);

% Tau = I* alpha
tau_aero = 1e-4; % Nm
omega_limit_rpm = 4000;
omega_limit = omega_limit_rpm * (2*pi/60); % rad/s

alpha_aero = tau_aero / mean_moi; % rad/s^2

% omega = alpha * time
time_to_limit = omega_limit / alpha_aero; % seconds

H = mean_moi * omega_limit; % N·m·s

fprintf('\nAngular acceleration due to 1e-4 Nm torque: %.4e rad/s^2\n', alpha_aero);
fprintf('Time to reach 4000 rpm: %.2f seconds\n', time_to_limit);
fprintf('Angular momentum capacity: %.3f N·m·s\n', H);