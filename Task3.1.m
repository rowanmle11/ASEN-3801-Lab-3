clear
clc
close all

%% bringing in files

filename = '2025_09_23_002_F10A10.csv.xlsx';
header_lines = 2;
data = readmatrix(filename, "NumHeaderLines", header_lines);

%% Angular Rate vs Time

input_data_time = data(:,1);
input_data_gyro_output = data(:,2);
input_data_rate_rpm = data(:,3);

data_time = input_data_time - input_data_time(1);
data_rate_rads = input_data_rate_rpm .* ((2.*pi)./60);

[P, S] = polyfit(data_rate_rads, input_data_gyro_output, 1);
K = P(1); % slope = adjusted scale factor
b = P(2); % intercept = bias

calib_rate = (input_data_gyro_output - b) ./ K;
f = polyval(P, data_rate_rads);

test_label = 'Freq = 1.0 Hz, Amp = 1.0 A';

K_mean = mean(K);
b_mean = mean(b);

figure
plot(data_time, data_rate_rads, '-r')
hold on
plot(data_time, input_data_gyro_output, '-b')
xlabel('Time (s)')
ylabel('Angular Rate (rad/s)')
legend('Angular Rate', 'Gyro Output')
title('Angular Rate vs Time')
grid on

%% Using all data files

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

    gyro_calib = (gyro - b)./K;

    [P, S] = polyfit(rate_rads, gyro, 1);
    K_all(i) = P(1);
    b_all(i) = P(2);

    calib_rate = (input_data_gyro_output - b_all) ./ K_all;

    figure
    scatter(rate_rads, gyro, 2, 'b', 'filled')
    hold on
    f = polyval(P, rate_rads);
    plot(rate_rads, f, 'r', 'LineWidth', 1.5)
    yline(b_mean, '--k')
    xlabel('Encoder Angular Rate (rad/s)')
    ylabel('Gyro Output (rad/s)')
    title(['Calibration: ', test_labels{i}])
    legend('Data', 'Linear Fit', 'Mean Bias', 'Location', 'best')
    grid on

    figure
    plot(time, rate_rads, 'b')
    hold on
    plot(time, gyro_calib, '--r')
    xlabel('Time (s)')
    ylabel('Angular Rate (rad/s)')
    legend('Encoder (Truth)', 'Calibrated Gyro', 'Location', 'best')
    title('Time History of Encoder vs Calibrated Gyro')
    grid on
end

K_mean = mean(K_all);
K_std = std(K_all);
b_mean = mean(b_all);
b_std = std(b_all);

%% 3.1c Part 3

T = table(test_labels', K_all', b_all', ...
           'VariableNames', {'Trial', 'ScaleFactor K', 'Bias b'});
disp(T)
fprintf('\nMean K = %.4f ± %.4f\n', K_mean, K_std)
fprintf('Mean b = %.4f ± %.4f\n', b_mean, b_std)

%% 3.1c Part 4

header_lines = 2;
filenames = {'2025_09_23_002_F02A05.csv.xlsx',...
             '2025_09_23_002_F04A10.csv.xlsx'};

K_vals = [-0.8712, -0.87181];
b_vals = [-0.0025171, -0.0038706];

for i = 1:numel(filenames)
    fname = filenames{i};
    K = K_vals(i);
    b = b_vals(i);

    data = readmatrix(fname, "NumHeaderLines", header_lines);
    time = data(:,1) - data(1,1);
    gyro_raw = data(:,2);
    rate_rpm = data(:,3);
    rate_true = rate_rpm .* (2.*pi)./60;

    calib_rate = (gyro_raw - b) ./ K;

    encoder_pos = cumtrapz(time, rate_true);
    gyro_pos = cumtrapz(time, calib_rate);

    rate_error = calib_rate - rate_true;
    pos_error = gyro_pos - encoder_pos;

    figure('Name', sprintf('Gyro Calibration Comparison: %s', fname), 'Units', 'normalized', 'Position', [0.1 0.1 0.75 0.8])
    subplot(3,2,1)
    plot(time, rate_true, 'b')
    hold on
    plot(time, calib_rate, '--r')
    xlabel('Time (s)')
    ylabel('Angular Rate (rad/s)')
    legend('Encoder (Truth)', 'Calibrated Gyro', 'Location', 'best')
    title('i: Angular Rate Comparison')
    grid on

    subplot(3,2,2)
    plot(time, rate_error, 'k')
    xlabel('Time (s)')
    ylabel('Rate Error (rad/s)')
    title('ii: Angular Rate Error')
    grid on

    subplot(3,2,3)
    plot(time, encoder_pos, 'b')
    plot(time, gyro_pos, '--r')
    xlabel('Time (s)')
    ylabel('Angular Position (rad)')
    legend('Encoder (Truth)', 'Gyro (Measured)', 'Location', 'best')
    title('iii: Angular Position Comparison')
    grid on

    subplot(3,2,4)
    plot(time, pos_error, 'm')
    xlabel('Time (s)')
    ylabel('Position Error (rad)')
    title('iv: Position Error vs Time')
    grid on

    subplot(3,2,[5,6])
    plot(rate_true, pos_error, '.')
    xlabel('Encoder Rate (rad/s)')
    ylabel('Position Error (rad)')
    title('v: Position Error vs Encoder Rate')
    grid on
end
