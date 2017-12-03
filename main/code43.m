% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2017, Andreas L. Fl?ten

%% Initialization and model definition
code42; % NB: Change this to the init file corresponding to your helicopter

global N mx

% Number of states and inputs
mx = size(A,2); % Number of states (number of columns in A)
mu = size(B,2); % Number of inputs(number of columns in B)

% Initial values
x0_lambda   = pi;                     
x0_r        = 0;                               
x0_p        = 0;                          
x0_p_dot    = 0 ;                     
x0_e        = 0;                
x0_e_dot    = 0;
x0 = [x0_lambda x0_r x0_p x0_p_dot x0_e x0_e_dot]'; 

% Time horizon and initialization
N  = 40;                                % Time horizon task 4
M  = N;                                 % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);                % Initialize z for the whole horizon
z0 = z;                                 % Initial value for optimization

% Generate the matrix Q and the vector c (objecitve function weights in the QP problem) 
Q1 = zeros(mx,mx);
Q1(1,1) = 1;                            % Weight on state x1
Q1(2,2) = 0;                            % Weight on state x2
Q1(3,3) = 0;                            % Weight on state x3
Q1(4,4) = 0;                            % Weight on state x4
Q1(5,5) = 0;                            % Weight on state x5
Q1(6,6) = 0;                            % Weight on state x6
q1 = 1;
q2 = 1;
P1 = hints_blkdiag(q1,q2);              % Weight on input
Q = 2*hints_genq2(Q1,P1,N,M,mu);        % Generate Q

%% Generate system matrixes for linear model
Aeq = hints_gena2(A,B,N,mx,mu);           % Generate A, hint: gena2
beq = [A*x0; zeros((N-1)*mx,1)];  


%Bounds
lowerBounds = [
    repmat([-Inf; -Inf; -Inf; -Inf; -Inf; -Inf], N, 1);
    repmat([-pi/6; -pi/3], N, 1)];

upperBounds = [ 
    repmat([Inf; Inf; Inf; Inf; Inf; Inf], N, 1);
    repmat([pi/6; pi/3], N, 1)];

%% Solve QP problem with linear model
% objective function
myf = @(lambda, p, e) ones(1, N)*(lambda.^2 + q1.*p.^2 + q2.*e.^2);
f = @(x) myf(x(1: 6: 6*N), x(3 : 6: 6*N), x(5: 6: 6*N));

options = optimset('Algorithm','active-set');
options = optimset(options, 'MaxFunEvals', 60000);
tic  
[z, fval] = fmincon(f,z0,[],[],Aeq,beq,lowerBounds,upperBounds,@func_constraint, options);
t1=toc;

%% Extract control inputs and states
%u1  = [z(N*mx+1:mu:N*mx+M*mu);z(N*mx+M*mu)]; % Control input from solution
%u2  = [z(N*mx+2:mu:N*mx+M*mu);z(N*mx+M*mu)];

u1 = z(N*6 + 1:2:N*8);
u2 = z(N*6 + 2:2:N*8);

%x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
%x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
%x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
%x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
%x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
%x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution

x1 = z(1:6:N*6);
x2 = z(2:6:N*6);
x3 = z(3:6:N*6);
x4 = z(4:6:N*6);
x5 = z(5:6:N*6);
x6 = z(6:6:N*6);

num_variables = 5/Theta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1  = [zero_padding; u1; zero_padding];
u2  = [zero_padding; u2; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; zero_padding];
x6  = [zero_padding; x6; zero_padding];

%% Plotting
t = 0:Theta_t:Theta_t*(length(u1)-1);

figure(43)
subplot(811)
stairs(t,u1),grid
ylabel('u1')
subplot(812)
stairs(t,u2),grid
ylabel('u2')
subplot(813)
plot(t,x1,'m',t,x1,'r'),grid
ylabel('lambda')
subplot(814)
plot(t,x2,'m',t,x2','r'),grid
ylabel('r')
subplot(815)
plot(t,x3,'m',t,x3,'r'),grid
ylabel('p')
subplot(816)
plot(t,x4,'m',t,x4','r'),grid
ylabel('pdot')
subplot(817)
plot(t,x5,'m',t,x5','r'),grid
ylabel('e')
subplot(818)
plot(t,x6,'m',t,x6','r'),grid
ylabel('edot')
hold off

input = [t' u1 u2];
x_star = [t' x1 x2 x3 x4 x5 x6];

ut1 = [t' u1];
ut2 = [t' u2];  
