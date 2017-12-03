% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2017, Andreas L. Fl?ten

%% Initialization and model definition
code31; % NB: Change this to the init file corresponding to your helicopter

%% Plotting

load('32.mat');

figure(32)
subplot(511)
stairs(t,u),grid
ylabel('u')
subplot(512)
hold on
lambda_prediction = plot(x_state(1,:),x_state(2,:),'b');
lambda_meassured = plot(x_state(1,:),x_state(6,:),'r');
legend([lambda_prediction, lambda_meassured], 'Lambda (predicted)', 'Lambda (meassured)');
hold off
ylabel('lambda') 
subplot(513)
hold on
lambda_prediction = plot(x_state(1,:),x_state(3,:),'b');
lambda_meassured = plot(x_state(1,:),x_state(7,:),'r');
legend([lambda_prediction, lambda_meassured], 'r (predicted)', 'r (meassured)');
hold off
ylabel('r')
subplot(514)
hold on
lambda_prediction = plot(x_state(1,:),x_state(4,:),'b');
lambda_meassured = plot(x_state(1,:),x_state(8,:),'r');
legend([lambda_prediction, lambda_meassured], 'p (predicted)', 'p (meassured)');
hold off
ylabel('p')
subplot(515)
hold on
lambda_prediction = plot(x_state(1,:),x_state(5,:),'b');
lambda_meassured = plot(x_state(1,:),x_state(9,:),'r');
legend([lambda_prediction, lambda_meassured], 'pdot (predicted)', 'pdot (meassured)');
hold off
xlabel('tid (s)'),ylabel('pdot')

%% For simulink
x_star = [t', x1, x2, x3, x4];
tu = [t' u]; 

