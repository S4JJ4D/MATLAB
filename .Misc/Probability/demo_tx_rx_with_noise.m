%% Example 8.10 of Bertsekas

close all;
clear;

% transmitted messages are of two types: a, b

R = 3;

a = R * [2, 1]/norm([2 1]);
b = R * [-3, 1]/norm([-3, 1]);

msg_space = [a;b];

quiver(0, 0, a(1), a(2), 'AutoScale', 'off', 'LineWidth', 1);
hold on;
quiver(0, 0, b(1), b(2), 'AutoScale', 'off', 'LineWidth', 1);

% draw the circle
theta_ = 0:.01:2*pi;
plot(R*cos(theta_), R*sin(theta_), 'r-');

axis([-5.5351    5.6028   -4.2339    4.5507]);

received_msg_q = quiver(0,0,nan,nan, 'AutoScale', 'off', 'LineWidth', 1);
received_msg_q.Color = 'k';

N = 100;

correct_estimate_count = 0;

tobj = title('Ratio of Correct Estimates = 0.0');

for i=1:N
    % transmitter sends a message
    msg_id = randi(2);
    transmitted_msg = msg_space(msg_id,:);

    % message is corrupted with noise during transit
    received_msg = transmitted_msg + randn(1, 2);
    set(received_msg_q, 'UData', received_msg(1), 'VData', received_msg(2));

    if dot(received_msg, a) - dot(received_msg, b) >=0 && all(transmitted_msg == a)
        correct_estimate_count = correct_estimate_count + 1;
    elseif dot(received_msg, b) - dot(received_msg, a) >=0 && all(transmitted_msg == b)
        correct_estimate_count = correct_estimate_count + 1;
    end

    tobj.String = ['Ratio of Correct Estimates = ', num2str(correct_estimate_count/i)];
    pause(.1);
end
