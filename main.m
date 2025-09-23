clear
clc
close all

filename = '2025_09_23_002_F02A05.csv.xlsx';
header_lines = 2;
data = readmatrix(filename, "NumHeaderLines", header_lines);

input_data_time = data(:,1);
input_data_gyro_output = data(:,2);
input_data_rate_rpm = data(:,3);

data_time = input_data_time - input_data_time(1);
data_rate_rads = input_data_rate_rpm .* ((2.*pi)./60);

figure
plot(data_time, data_rate_rads, '-r')
hold on
plot(data_time, input_data_gyro_output, '-b')
hold on
xlabel('Time (s)')
ylabel('Angular Rate (rad/s)')
legend('Angular Rate', 'Gyro Output')
grid on

figure
plot(data_rate_rads, input_data_gyro_output)
hold on
title('Gyro Output vs Angular Rate')
xlabel('Angular Rate (rad/s)')
ylabel('Gyro Output (rad/s)')
grid on

[P, S] = polyfit(data_rate_rads, input_data_gyro_output, 1);
f = polyval(P, data_rate_rads);

plot(f, P) % incomplete

%% For shits and giggles

fprintf('Best Regards,\n')
fprintf('Lab 3 Group 4')