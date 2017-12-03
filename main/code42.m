% Initialization for the helicopter assignment in TTK4135.
% Run this file before you execute QuaRC -> Build.

% Updated spring 2017, Andreas L. Flåten

init03;

%Code task 4.1
I = eye(6);
A_c = [0 1 0 0 0 0; 0 0 -K_2 0 0 0; 0 0 0 1 0 0; 0 0 -K_1*K_pp -K_1*K_pd 0 0; 0 0 0 0 0 1; 0 0 0 0 -K_3*K_ep -K_3*K_ed];
B_c = [0 0; 0 0; 0 0; K_1*K_pp 0; 0 0; 0 K_3*K_ep];

%Code task 4.2
Theta_t = 0.25;
A = I + Theta_t*A_c;
B = Theta_t*B_c;

