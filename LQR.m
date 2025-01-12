%%  Clear Workspace
clc; clear; close all

%% Parameters

[m , M, h, a, r, I, I_x, I_y, I_z, g] = deal(0.141, 2.313, 0.254, 0.165, 0.0615, 2.5e-4, 0.1986, 0.1942, 0.0056, 9.81);

%% Linear State Space Realisation
A = [0, 0, 0, 1, 0, 0;...
       0, 0, 0, 0, 1, 0;...
       0, 0, 0, 0, 0, 1;...
       0, 0, -1*(M^2*g*h^2*r^2)/(2*I*I_y + 2*I*M*h^2 + I_y*M*r^2 + 2*I_y*m*r^2 + 2*M*h^2*m*r^2), 0, 0, 0;...
       0, 0, 0, 0, 0, 0;...
       0, 0, (2*g*h*M^2*r^2 + 4*g*h*m*M*r^2 + 4*I*g*h*M)/(4*I*I_y + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 + 4*M*h^2*m*r^2), 0, 0, 0];

B = [0, 0;...
      0, 0;...
      0, 0;...
      (r*(M*h^2 + M*r*h + I_y))/(2*I*I_y + 2*I*M*h^2 + I_y*M*r^2 + 2*I_y*m*r^2 + 2*M*h^2*m*r^2), (r*(M*h^2 + M*r*h + I_y))/(2*I*I_y + 2*I*M*h^2 + I_y*M*r^2 + 2*I_y*m*r^2 + 2*M*h^2*m*r^2);...
    -1*(a*r)/(2*m*a^2*r^2 + 2*I*a^2 + I_z*r^2), (a*r)/(2*m*a^2*r^2 + 2*I*a^2 + I_z*r^2);...
    -1*(4*I + 2*M*r^2 + 4*m*r^2 + 2*M*h*r)/(4*I*I_y + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 + 4*M*h^2*m*r^2), -1*(4*I + 2*M*r^2 + 4*m*r^2 + 2*M*h*r)/(4*I*I_y + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 + 4*M*h^2*m*r^2)];

C = eye(6);

D = zeros(6, 2);

states = {'rho' 'psi' 'theta' 'rho_dot' 'psi_dot' 'theta_dot'};
inputs = {'T_l' 'T_r'};
outputs = {'rho' 'psi' 'theta', 'rho_dot' 'psi_dot' 'theta_dot'};

sys_ol = ss(A, B, C, D, 'statename', states, 'inputname', inputs,'outputname', outputs);
G = tf(sys_ol);
display(G)

%% Compute Pole Placing Gain
% Q chosen for faster response in rho by higher penalization
% and a gradual transition to linear motion by penalizing psi less
Q = diag([10, 0.01, 15, 1, 0.1, 1]); 
% R chosen to ensure 
R = 0.001*eye(2);
K = lqr(A, B, Q, R, 0);
display(K)

%% Check Placed Eigenvalues
poles = eig(A - B*K);
display(poles)

%% Closed Loop
sys_cl = ss(A - B*K, B, C - D*K, D, 'statename', states, 'inputname', inputs, 'outputname', outputs);
display(sys_cl)
%% Torque Inputs

simulationTime = 20; % seconds
steps = 10000;
timeHorizon = linspace(0, simulationTime, steps);
T_l = 0 * (heaviside(timeHorizon) - heaviside(timeHorizon - 0.1));
T_r = 0 * (heaviside(timeHorizon) - heaviside(timeHorizon - 0.11));
u = [T_l;...
      T_r];


%% Linear Simulation

initialCondition = [1; pi/4; pi/6; zeros(3, 1)];
[y, tOut] = lsim(sys_cl, u, timeHorizon, initialCondition);

%% Nonlinear Simulation settings

Solver = "ODE45";
Options = odeset("RelTol", 1e-6, "AbsTol", 1e-6);

%% Nonlinear Simulation

[~, X] = ode45(@(t, X) RHS(t, X, timeHorizon, K, m, M, h, a, r, I, I_x, I_y , I_z, g), timeHorizon, initialCondition, Options);

%% Visualize Results

figure
sgtitle(" LQR and Nonlinear Closed Loop Simulation Results")
subplot(6, 1, 1)
plot(timeHorizon, X(:, 1), "DisplayName", "Nonlinear Feedback")
hold on
plot(tOut, y(:, 1), "DisplayName", "LQR")
legend("location", "best")
title("Radial Distance ($\rho$)", "Interpreter", "latex")

subplot(6, 1, 2)
plot(timeHorizon, X(:, 2))
hold on
plot(tOut, y(:, 2))
title("Angular Offset ($\psi$)", "Interpreter", "latex")

subplot(6, 1, 3)
plot(timeHorizon, X(:, 3))
hold on
plot(tOut, y(:, 3))
title("Pendulum Tilt Angle ($\theta$)", "Interpreter", "latex")

subplot(6, 1, 4)
plot(timeHorizon, X(:, 4))
hold on
plot(tOut, y(:, 4))
title("Linear Acceleration ($\dot\rho$)", "Interpreter", "latex")

subplot(6, 1, 5)
plot(timeHorizon, X(:, 5))
hold on
plot(tOut, y(:, 5))
title("Angular Acceleration ($\dot\psi$)", "Interpreter", "latex")

subplot(6, 1, 6)
plot(timeHorizon, X(:, 6))
hold on
plot(tOut, y(:, 6))
title("Pendulum Angular Acceleration ($\dot\theta$)", "Interpreter", "latex")
xlabel("$t$", "Interpreter", "latex")

%% Segway Path in Cartesian  Coordinates
% plot trajectories
x_nonlinear = X(:, 1) .* cos(X(:, 2));
y_nonlinear = X(:, 1) .* sin(X(:, 2));

x_linear = y(:, 1) .* cos(y(:, 2));
y_linear = y(:, 1) .* sin(y(:, 2));

figure
hold on
grid on
title("Segway Path in Cartesian Coordinates")
xlabel("$x(t)$", "Interpreter", "latex")
ylabel("$y(t)$", "Interpreter", "latex")

% animation
plot(x_nonlinear, y_nonlinear, "--", "DisplayName", "Nonlinear Path")
plot(x_linear, y_linear, "--", "DisplayName", "LQR Path")
legend("location", "best")

% add markers for linear and nonlinear trajectories
nonlinearMarker = plot(0, 0, "bo", "MarkerSize", 6, "MarkerFaceColor", "b","DisplayName", "Nonlinear Segway");
LQRMarker = plot(0, 0, "ro", "MarkerSize", 6, "MarkerFaceColor", "r","DisplayName", "LQR Segway");
animationPauseTime = simulationTime / steps;
tic
for i = 1:8:length(timeHorizon)
    % update markers for linear and nonlinear positions
    LQRMarker.XData = x_linear(i);
    LQRMarker.YData = y_linear(i);
    
    nonlinearMarker.XData = x_nonlinear(i);
    nonlinearMarker.YData = y_nonlinear(i);
    
    % pause for animation speed control
    pause(animationPauseTime);
end
animationTime = toc;
fprintf(" animationTime = %.2f seconds", animationTime)
%% ODE Function

function dxdt = RHS(~, X, ~, K, m, M, h, a, r, I, I_x, I_y, I_z, g)
    % state variables
    rho = X(1);
    psi = X(2);
    theta = X(3);
    rho_dot = X(4);
    psi_dot = X(5);
    theta_dot = X(6);

    % state vector
    x = [rho; psi; theta; rho_dot; psi_dot; theta_dot];

    % control input (state feedback law)
    u = -K * x; 
    T_l = u(1); % left wheel torque
    T_r = u(2); % right wheel torque

    % dynamical state equations
    dxdt = zeros(6, 1); % preallocation for efficiency
    dxdt(1) = rho_dot;
    dxdt(2) = psi_dot;
    dxdt(3) = theta_dot;
    dxdt(4) = (r * (I_y * T_l + I_y * T_r + M * T_l * h^2 + M * T_r * h^2 + M^2 * h^3 * r * theta_dot^2 * sin(theta) ...
        + M * T_l * h * r * cos(theta) + M * T_r * h * r * cos(theta) + M^2 * h^3 * psi_dot^2 * r * sin(theta) ...
        - M^2 * h^3 * r * rho_dot^2 * cos(theta)^2 * sin(theta) - M^2 * g * h^2 * r * cos(theta) * sin(theta) ...
        + I_y * M * h * psi_dot^2 * r * sin(theta) + I_y * M * h * r * theta_dot^2 * sin(theta) ...
        - I_x * M * h * r * rho_dot^2 * cos(theta)^2 * sin(theta) + I_z * M * h * r * rho_dot^2 * cos(theta)^2 * sin(theta))) ...
        / (-M^2 * h^2 * r^2 * cos(theta)^2 + M^2 * h^2 * r^2 + 2 * m * M * h^2 * r^2 + 2 * I * M * h^2 + I_y * M * r^2 + 2 * I_y * m * r^2 + 2 * I * I_y);
    dxdt(5) = -1 * (r * (m * psi_dot * r * theta_dot * sin(2 * theta) * h^2 + M * psi_dot * r * rho_dot * sin(theta) * h ...
        + T_l * a - T_r * a + I_x * psi_dot * r * theta_dot * sin(2 * theta) - I_z * psi_dot * r * theta_dot * sin(2 * theta))) ...
        / (2 * m * a^2 * r^2 + 2 * I * a^2 + M * h^2 * r^2 * sin(theta)^2 + I_z * r^2 * cos(theta)^2 + I_x * r^2 * sin(theta)^2);
    dxdt(6) = -1 * (4 * I * T_l + 4 * I * T_r + 2 * M * T_l * r^2 + 2 * M * T_r * r^2 + 4 * T_l * m * r^2 + 4 * T_r * m * r^2 ...
        - 2 * I * I_x * rho_dot^2 * sin(2 * theta) + 2 * I * I_z * rho_dot^2 * sin(2 * theta) - 2 * I_x * m * r^2 * rho_dot^2 * sin(2 * theta) ...
        + 2 * I_z * m * r^2 * rho_dot^2 * sin(2 * theta) - 2 * M^2 * g * h * r^2 * sin(theta) - 4 * I * M * g * h * sin(theta) ...
        + 2 * M * T_l * h * r * cos(theta) + 2 * M * T_r * h * r * cos(theta) + M^2 * h^2 * psi_dot^2 * r^2 * sin(2 * theta) ...
        - M^2 * h^2 * r^2 * rho_dot^2 * sin(2 * theta) + M^2 * h^2 * r^2 * theta_dot^2 * sin(2 * theta) - 2 * I * M * h^2 * rho_dot^2 * sin(2 * theta) ...
        - I_x * M * r^2 * rho_dot^2 * sin(2 * theta) + I_z * M * r^2 * rho_dot^2 * sin(2 * theta) - 4 * M * g * h * m * r^2 * sin(theta) ...
        - 2 * M * h^2 * m * r^2 * rho_dot^2 * sin(2 * theta)) ...
        / (4 * I * I_y + M^2 * h^2 * r^2 + 4 * I * M * h^2 + 2 * I_y * M * r^2 + 4 * I_y * m * r^2 + 4 * M * h^2 * m * r^2 - M^2 * h^2 * r^2 * cos(2 * theta));
end
