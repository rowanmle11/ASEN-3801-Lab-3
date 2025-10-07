% Lab 3.3 Spacecraft MOI
clc;
clear;
close all;


% Lab 3.2 Code

trial_torque4 = "2025_09_30_002_SC_T4t5.xlsx";
trial_torque8 = "2025_09_30_002_SC_T8t5.xlsx";
trial_torque12 = "2025_09_30_002_SC_T12t5.xlsx";
trial_torque16 = "2025_09_30_002_SC_T16t5.xlsx";
trial_torque20 = "2025_09_30_002_SC_T20t5.xlsx";

trials = [trial_torque4 trial_torque8 trial_torque12 trial_torque16 trial_torque20];
trials_data = cell(1, length(trials));
torques = [4 8 12 16 20];

trials = [trial_torque4 trial_torque8 trial_torque12 trial_torque16 trial_torque20];
trials_data = cell(1, length(trials));
torques = [4 8 12 16 20];

figure(1);
sgtitle('Torque vs Time for All Trials');

figure(2);
sgtitle('Magnitude of Angular Velocity vs Time for All Trials')

for i = 1:length(trials)

    trials_data{i} = readmatrix(trials(i));
    time{i} = (trials_data{i}(:,1))/(10^3); %conversion from ms to sec
    torque{i} = (trials_data{i}(:,4))*(25.5); %multiply by torque constant mNm
    angular_velocity{i} = (-1*trials_data{i}(:,3))* ((2*pi)/60); %Convert from RPM to rad/s

    t_min = 1;
    t_max = 6;
    
    %Polyfit for alpha
    linear_idx = (time{i} >= t_min & time{i} <= t_max);
    p = polyfit(time{i}(linear_idx), angular_velocity{i}(linear_idx), 1);
    w_fit = polyval(p, time{i}(linear_idx));

    %Polyfit for Torque
    linear_idx_tor = (time{i} >= t_min & time{i} <= t_max);
    avg_torque = mean(torque{i}(linear_idx_tor));

    figure(1)
    subplot(2,3,i)
    plot(time{i}, torque{i})
    hold on;
    yline(avg_torque,'r', 'LineWidth', 2)
    title(['Torque for Trial ' num2str(i) ': ' num2str(torques(i)) ' mNm']);
    xlabel('Time (s)');
    ylabel('Torque (mNm)');
    legend('Measured Torque', 'Average Torque', 'Location', 'southwest');
    grid on;
    hold off; % Release the hold on the current figure

    figure(2)
    subplot(2,3,i)
    plot(time{i}, angular_velocity{i}, 'b');
    hold on;
    plot(time{i}(linear_idx), w_fit, 'r', 'LineWidth', 2);
    title(['Angular Velocity for Trial ' num2str(i) ': ' num2str(torques(i)) ' Nm']);
    xlabel('Time (s)');
    ylabel('Angular Velocity (rad/s)');
    legend('Measured Data', 'Linear Fit', 'Location', 'best');
    grid on;

    avg_torque(i) = avg_torque;

    alpha(i) = p(1);
    moi(i) = (avg_torque(i)*(0.001))./ alpha(i); % Using average torque and angular acceleration
    fprintf('Trial %d: Torque = %.3f mNm | Alpha = %.3f rad/s^2 | I = %.4e kg*m^2\n', ...
        i, avg_torque(i), alpha(i), moi(i));
end
mean_alpha = mean(alpha(i));
stdev = std(moi);
mean_moi = mean(moi); 

fprintf('Standard Deviation of Moment of Inertia: %.4e kg*m^2\n', stdev);
fprintf('Mean Moment of Inertia: %.4e kg*m^2\n', mean_moi);