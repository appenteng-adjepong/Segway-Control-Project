%%  Clear Workspace
clc; clear; close all

%% Parameters

[m , M, h, a, r, I, I_x, I_y, I_z, g] = deal(0.141, 2.313, 0.254, 0.165, 0.0615, 2.5e-4, 0.1986, 0.1942, 0.0056, 9.81);

%% Torque Inputs

simulationTime = 20; % seconds
steps = 10000;
timeHorizon = linspace(0, simulationTime, steps);
T_l = 1 * (heaviside(timeHorizon) - heaviside(timeHorizon - 0.1));
T_r = 1 * (heaviside(timeHorizon) - heaviside(timeHorizon - 0.11));

% visualize inputs
figure
sgtitle("Torque Inputs on Left and Right Wheels")
subplot(2, 1, 1)
plot(timeHorizon, T_l, "b", "LineWidth", 1)
xlim([0, 1])
xlabel("$t$", "Interpreter", "latex")
ylabel ("$T_l(t)$", "Interpreter", "latex")
subplot(2, 1, 2)
plot(timeHorizon, T_r, "r", "LineWidth", 1)
xlim([0, 1])
xlabel("$t$", "Interpreter", "latex")
ylabel("$T_r(t)$", "Interpreter", "latex")

%% Simulation settings

Solver = "ODE45";
Options = odeset("RelTol", 1e-6, "AbsTol", 1e-6);

%% Open Loop Nonlinear Simulation

initialCondition = [zeros(1, 6)];
[~, X] = ode45(@(t, X) RHS(t, X, timeHorizon, T_l, T_r, m, M, h, a, r, I, I_x, I_y , I_z, g), timeHorizon, initialCondition, Options);

%% Visualize Results

figure
sgtitle("Nonlinear Simulation Results")
subplot(6, 1, 1)
plot(timeHorizon, X(:, 1))
title("Radial Distance ($\rho$)", "Interpreter", "latex")

subplot(6, 1, 2)
plot(timeHorizon, X(:, 2))
title("Angular Offset ($\psi$)", "Interpreter", "latex")

subplot(6, 1, 3)
plot(timeHorizon, X(:, 3))
title("Pendulum Tilt Angle ($\theta$)", "Interpreter", "latex")

subplot(6, 1, 4)
plot(timeHorizon, X(:, 4))
title("Linear Acceleration ($\dot\rho$)", "Interpreter", "latex")

subplot(6, 1, 5)
plot(timeHorizon, X(:, 5))
title("Angular Acceleration ($\dot\psi$)", "Interpreter", "latex")

subplot(6, 1, 6)
plot(timeHorizon, X(:, 6))
title("Pendulum Angular Acceleration ($\dot\theta$)", "Interpreter", "latex")
xlabel("$t$", "Interpreter", "latex")

%% Segway Path in Cartesian  Coordinates

x = X(:, 1) .* cos(X(:, 2));
y = X(:, 1) .* sin(X(:, 2));

figure
hold on
grid on
plot(x, y, "--", "DisplayName", "Path")
title("Segway Path in Cartesian Coordinates")
xlabel("$x$", "Interpreter", "latex")
ylabel("$y$", "Interpreter", "latex")

% animation
marker = plot(0, 0, "bo", "MarkerSize", 6, "MarkerFaceColor", "b","DisplayName", "Segway");
legend(Location="best")
animationPauseTime = simulationTime / steps;
tic
for i = 1:12:length(timeHorizon)
    % update markers for segway positions
    
    marker.XData = x(i);
    marker.YData = y(i);
    
    % pause for animation speed control
    pause(animationPauseTime);
end
animationTime = toc;
fprintf(" animationTime = %.2f seconds", animationTime)
%% ODE Function
function dxdt = RHS(t, X, timeHorizon, T_l, T_r, m, M, h, a, r, I, I_x, I_y , I_z, g)
    T_l = interp1(timeHorizon, T_l, t);
    T_r = interp1(timeHorizon, T_r , t);

   % state variables
    theta = X(3);
    rho_dot = X(4);
    psi_dot = X(5);
    theta_dot = X(6);

% dynamical state equations
    dxdt = zeros(6, 1); % preallocation for efficiency
    dxdt(1) = rho_dot;
    dxdt(2) = psi_dot;
    dxdt(3) = theta_dot;
    dxdt(4) = (r*(I_y*T_l + I_y*T_r + M*T_l*h^2 + M*T_r*h^2 + M^2*h^3*r*theta_dot^2*sin(theta) ...
         + M*T_l*h*r*cos(theta) + M*T_r*h*r*cos(theta) + M^2*h^3*psi_dot^2*r*sin(theta) ...
         - M^2*h^3*r*rho_dot^2*cos(theta)^2*sin(theta) - M^2*g*h^2*r*cos(theta)*sin(theta) ...
         + I_y*M*h*psi_dot^2*r*sin(theta) + I_y*M*h*r*theta_dot^2*sin(theta) ...
         - I_x*M*h*r*rho_dot^2*cos(theta)^2*sin(theta) + I_z*M*h*r*rho_dot^2*cos(theta)^2*sin(theta))) ...
         / (- M^2*h^2*r^2*cos(theta)^2 + M^2*h^2*r^2 + 2*m*M*h^2*r^2 + 2*I*M*h^2 + I_y*M*r^2 ...
         + 2*I_y*m*r^2 + 2*I*I_y);
    dxdt(5 )= -(r*(m*psi_dot*r*theta_dot*sin(2*theta)*h^2 + M*psi_dot*r*rho_dot*sin(theta)*h + T_l*a ...
         - T_r*a + I_x*psi_dot*r*theta_dot*sin(2*theta) - I_z*psi_dot*r*theta_dot*sin(2*theta))) ...
         / (2*m*a^2*r^2 + 2*I*a^2 + M*h^2*r^2*sin(theta)^2 + I_z*r^2*cos(theta)^2 + I_x*r^2*sin(theta)^2);
    dxdt(6) = -(4*I*T_l + 4*I*T_r + 2*M*T_l*r^2 + 2*M*T_r*r^2 + 4*T_l*m*r^2 + 4*T_r*m*r^2 ...
         - 2*I*I_x*rho_dot^2*sin(2*theta) + 2*I*I_z*rho_dot^2*sin(2*theta) ...
         - 2*I_x*m*r^2*rho_dot^2*sin(2*theta) + 2*I_z*m*r^2*rho_dot^2*sin(2*theta) ...
         - 2*M^2*g*h*r^2*sin(theta) - 4*I*M*g*h*sin(theta) + 2*M*T_l*h*r*cos(theta) ...
         + 2*M*T_r*h*r*cos(theta) + M^2*h^2*psi_dot^2*r^2*sin(2*theta) ...
         - M^2*h^2*r^2*rho_dot^2*sin(2*theta) + M^2*h^2*r^2*theta_dot^2*sin(2*theta) ...
         - 2*I*M*h^2*rho_dot^2*sin(2*theta) - I_x*M*r^2*rho_dot^2*sin(2*theta) ...
         + I_z*M*r^2*rho_dot^2*sin(2*theta) - 4*M*g*h*m*r^2*sin(theta) ...
         - 2*M*h^2*m*r^2*rho_dot^2*sin(2*theta)) ...
         / (4*I*I_y + M^2*h^2*r^2 + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 ...
         + 4*M*h^2*m*r^2 - M^2*h^2*r^2*cos(2*theta));
end