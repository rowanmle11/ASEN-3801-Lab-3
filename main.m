clear
clc
close all

%% bringing in files

filename = '2025_09_23_002_F02A05.csv.xlsx';
header_lines = 2;
data = readmatrix(filename, "NumHeaderLines", header_lines);

%% Angular Rate vs Time

input_data_time = data(:,1);
input_data_gyro_output = data(:,2);
input_data_rate_rpm = data(:,3);

data_time = input_data_time - input_data_time(1);
data_rate_rads = input_data_rate_rpm .* ((2.*pi)./60);

figure
plot(data_time, data_rate_rads, '-r')
hold on
plot(data_time, input_data_gyro_output, '-b')
xlabel('Time (s)')
ylabel('Angular Rate (rad/s)')
legend('Angular Rate', 'Gyro Output')
grid on

% [P, S] = polyfit(data_rate_rads, input_data_gyro_output, 1);
% K = P(1); % slope = adjusted scale factor
% b = P(2); % intercept = bias

% calib_rate = (input_data_gyro_output - b) ./ K;
% f = polyval(P, data_rate_rads);

files = {'2025_09_23_002_F02A05.csv.xlsx',...
         '2025_09_23_002_F04A10.csv.xlsx',...
         '2025_09_23_002_F10A10.csv.xlsx'};

test_labels = {'Freq = 0.2 Hz, Amp = 0.5 A',...
               'Freq = 0.4 Hz, Amp = 1.0 A',...
               'Freq = 1.0 Hz, Amp = 1.0 A'};

K_all = [];
b_all = [];

for i = 1:length(files)
    data = readmatrix(files{i}, "NumHeaderLines", header_lines);
    time = data(:,1) - data(1,1);
    gyro = data(:,2);
    rate_rads = data(:,3) .* (2.*pi./60);

    [P, S] = polyfit(rate_rads, gyro, 1);
    K_all(i) = P(1);
    b_all(i) = P(2);

    calib_rate = (input_data_gyro_output - b_all) ./ K_all;

    figure
    scatter(rate_rads, gyro, 2, 'b', 'filled')
    hold on
    f = polyval(P, rate_rads);
    plot(rate_rads, f, 'r', 'LineWidth', 1.5)
    xlabel('Encoder Angular Rate (rad/s)')
    ylabel('Gyro Output (rad/s)')
    title(['Calibration: ', test_labels{i}])
    legend('Data', 'Linear Fit', 'Location', 'best')
    grid on

end

K_mean = mean(K_all);
K_std = std(K_all);
b_mean = mean(b_all);
b_std = std(b_all);

%% 3.1d (Gyro Analysis)

dt = mean(diff(data_time));
encoder_pos = cumtrapz(data_time, data_rate_rads);
gyro_pos = cumtrapz(data_time, calib_rate);

rate_error = calib_rate - data_rate_rads;
pos_error = gyro_pos - encoder_pos;

%% PLOTS

figure
plot(data_time, data_rate_rads)
hold on
plot(data_time, calib_rate, '--r')
xlabel('Time (s)')
ylabel('Angular Rate (rad/s)')
legend('Encoder Rate (Truth)', 'Calibrated Gyro Rate')
title('Time History of Calibrated Gyro vs Encoder')
grid on

T = table(test_labels', K_all', b_all', ...
           'VariableNames', {'Trial', 'ScaleFactor K', 'Bias b'});
disp(T)
fprintf('\nMean K = %.4f ± %.4f\n', K_mean, K_std)
fprintf('Mean b = %.4f ± %.4f\n', b_mean, b_std)

figure
plot(data_time, rate_error)
xlabel('Time (s)')
ylabel('Rate Error (rad/s)')
title('Angular Rate Measurement Error')
legend('0.2 Hz, 0.5 A', '0.4 Hz, 1.0 A', '1.0 Hz, 1.0 A')
grid on

figure
plot(data_time, encoder_pos, data_time, gyro_pos)
xlabel('Time (s)')
ylabel('Angular Position (rad)')
legend('Encoder (Truth)', 'Gyro (Measured)')
title('Angular Position Comparison')
legend('0.2 Hz, 0.5 A', '0.4 Hz, 1.0 A', '1.0 Hz, 1.0 A')
grid on

figure
plot(data_time, pos_error)
xlabel('Time (s)')
ylabel('Angular Position Error (rad)')
title('Angular Position Error vs. Time')
legend('0.2 Hz, 0.5 A', '0.4 Hz, 1.0 A', '1.0 Hz, 1.0 A')
grid on

figure
plot(data_rate_rads, pos_error, '.')
xlabel('Encoder Rate (rad/s)')
ylabel('Angular Position Error (rad)')
title('Angular Position Error vs Encoder Rate')
legend('0.2 Hz, 0.5 A', '0.4 Hz, 1.0 A', '1.0 Hz, 1.0 A')
grid on
