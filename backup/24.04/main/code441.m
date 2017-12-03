% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2017, Andreas L. Fl?ten

%% Initialization and model definition
code43; % NB: Change this to the init file corresponding to your helicopter

%% Plotting
load('441.mat');

figure(441)
subplot(811)
stairs(t,u1),grid
ylabel('u1')
subplot(812)
stairs(t,u2),grid
ylabel('u2')
subplot(813)
hold on
lambda_prediction = plot(x_state(1,:),x_state(2,:),'b'), grid;
lambda_meassured = plot(x_state(1,:),x_state(8,:),'r');
legend([lambda_prediction, lambda_meassured], 'predicted', 'meassured');
hold off
ylabel('lambda')
subplot(814)
hold on
plot(x_state(1,:),x_state(3,:),'b'), grid;
plot(x_state(1,:),x_state(9,:),'r');
hold off
ylabel('r')
subplot(815)
hold on
plot(x_state(1,:),x_state(4,:),'b'), grid;
plot(x_state(1,:),x_state(10,:),'r');
hold off
ylabel('p')
subplot(816)
hold on
plot(x_state(1,:),x_state(5,:),'b'), grid;
plot(x_state(1,:),x_state(11,:),'r');
hold off
ylabel('pdot')
subplot(817)
hold on
plot(x_state(1,:),x_state(6,:),'b'), grid;
plot(x_state(1,:),x_state(12,:),'r');
hold off
ylabel('e')
subplot(818)
hold on
plot(x_state(1,:),x_state(7,:),'b'), grid;
plot(x_state(1,:),x_state(13,:),'r');
hold off
ylabel('edot')
xlabel('tid (s)') 

%%For simulink
input = [t' u1 u2];
x_star = [t' x1 x2 x3 x4 x5 x6];
